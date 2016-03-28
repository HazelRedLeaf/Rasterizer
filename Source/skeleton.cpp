#include <iostream>
#include <algorithm>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::vec2;
using glm::vec3;
using glm::mat3;
using glm::ivec2;



/* ----------------------------------------------------------------------------*/
/* STRUCTS                                                                     */

struct Pixel {
	int x;
	int y;
	float zinv;
	vec3 pos3d;
};

struct Vertex {
	vec3 position;
};


/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
int t;
vector<Triangle> triangles;
vec3 currentColor;
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];

// camera variables
vec3 cameraPos(0.f, 0.f, -3.001f);
float f = 500.f;
float cameraSpeed = 0.002f;
float yaw = 0; // Yaw angle controlling camera rotation around y-axis
mat3 R(vec3( 0, 0, 1),
       vec3( 0, 1, 0),
       vec3(-1, 0, 0));

// light variables
float lightSpeed = 0.2f;
vec3 lightPos(0,-0.5,-0.7);
vec3 lightPower = 11.1f*vec3(1, 1, 1);
vec3 indirectLightPowerPerArea = 0.5f*vec3(1, 1, 1);
vec3 currentNormal;
vec3 currentReflectance;


/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw();
void updateCameraAngle(float angle);
void VertexShader(const Vertex& v, Pixel& p);
void PixelShader(const Pixel& p);
void Interpolate(Pixel a, Pixel b, vector<Pixel>& result);
void ComputePolygonRows(const vector<Pixel>& vertexPixels, 
	vector<Pixel>& leftPixels,
	vector<Pixel>& rightPixels);
void DrawPolygonRows(const vector<Pixel>& leftPixels,
					 const vector<Pixel>& rightPixels);
void DrawPolygon(const vector<Vertex>& vertices);

int main(int argc, char* argv[])
{
    //is necessary for multithreaded access
    //XInitThreads();

	screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT);
	t = SDL_GetTicks();	// Set start value for timer.

	//Load scene triangles
    LoadTestModel(triangles);
    //initialize camera angle with default yaw
    updateCameraAngle(yaw);

	while(NoQuitMessageSDL())
	{
		Update();
		Draw();
	}

	SDL_SaveBMP( screen, "screenshot.bmp" );
	return 0;
}

void Update()
{
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2-t);
	t = t2;
	cout << "Render time: " << dt << " ms." << endl;

	Uint8* keystate = SDL_GetKeyState(0);
    //Move camera
    if(keystate[SDLK_UP]) {
        cameraPos.z += cameraSpeed;
    } 
    if(keystate[SDLK_DOWN]) {
        cameraPos.z -= cameraSpeed;
    } 
    if(keystate[SDLK_LEFT]) {
        cameraPos.x -= cameraSpeed;
        updateCameraAngle(-M_PI/18.f);
    } 
    if(keystate[SDLK_RIGHT]) {
        cameraPos.x += cameraSpeed;
        updateCameraAngle(M_PI/18.f);
    }
    //Move light source
    if(keystate[SDLK_w]) {
        lightPos.z += lightSpeed;
    }
    if(keystate[SDLK_s]){
        lightPos.z -= lightSpeed;
    }
    if(keystate[SDLK_d]){
        lightPos.x -= lightSpeed;
    }
    if(keystate[SDLK_a]){
        lightPos.x += lightSpeed;
    }
    if(keystate[SDLK_q]){
        lightPos.y += lightSpeed;
    }
    if(keystate[SDLK_e]){
        lightPos.y -= lightSpeed;
    }
}

void Draw()
{
	SDL_FillRect(screen, 0, 0);

	if(SDL_MUSTLOCK(screen))
		SDL_LockSurface(screen);

	for(int y=0; y<SCREEN_HEIGHT; ++y)
		for(int x=0; x<SCREEN_WIDTH; ++x)
			depthBuffer[y][x] = 0;

	for(size_t i = 0; i < triangles.size(); ++i)
	{
		vector<Vertex> vertices(3);
		vertices[0].position = triangles[i].v0;
		vertices[1].position = triangles[i].v1;
		vertices[2].position = triangles[i].v2;
		currentNormal = triangles[i].normal;
		currentReflectance = triangles[i].color;

		DrawPolygon(vertices);
	}

	if(SDL_MUSTLOCK(screen))
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}

void VertexShader(const Vertex& v, Pixel& p)
{
	// column vectors from rotation matrix
	vec3 R1 (R[0][0], R[1][0], R[2][0]);
	vec3 R2 (R[0][1], R[1][1], R[2][1]);
	vec3 R3 (R[0][2], R[1][2], R[2][2]);

	// get position of point in the camera's coordinate system
	float X = v.position.x - cameraPos.x;
	float Y = v.position.y - cameraPos.y;
	float Z = v.position.z - cameraPos.z;
	//vec3 P (X, Y, Z);
	X = X * R1.x + X * R1.z;
	Z = Z * R3.x + Z * R3.z;

	// store 1/Z for the depth buffer
	p.zinv = 1/Z;
	// project (X,Y,Z) to (x,y,f)
	p.x = int(f*X/Z) + SCREEN_WIDTH/2;
	p.y = int(f*Y/Z) + SCREEN_HEIGHT/2;	
	// store the 3D position of the vertex to the corresponding
	// variable in Pixel
	p.pos3d = v.position;
}

void PixelShader(const Pixel& p)
{
	int x = p.x;
	int y = p.y;

	// the current pixel is closer to the camera than what was drawn
	// at this position before (or infinity - nothing was drawn)
	if(p.zinv > depthBuffer[y][x]) {
		depthBuffer[y][x] = p.zinv;

		// vec3 D = lightPower*max(glm::dot(u_r,u_n),0)/4*M_PI*rsq;
		//Direction from surface point to light source
    	vec3 r = lightPos - p.pos3d;
    	//distance between intersection and light source
    	float dist = glm::length(r);
    	//distance squared between intersection and light source
    	float rsq = glm::dot(r, r);
    	vec3 u_r = glm::normalize(r);
    	vec3 u_n = currentNormal;
		float D_prime = max(glm::dot(u_r,u_n),0.f)/((float)4*M_PI*rsq);
		vec3 D = lightPower*D_prime;

		// R = rho * (D + N)
		vec3 R = currentReflectance * (D + indirectLightPowerPerArea);
		PutPixelSDL(screen, x, y, R);
	}
}

void Interpolate(Pixel a, Pixel b, vector<Pixel>& result)
{
	int N = result.size();
	float stepX = (b.x - a.x) / float(max(N-1,1));
	float stepY = (b.y - a.y) / float(max(N-1,1));
	float stepZ = (b.zinv - a.zinv) / float(max(N-1,1));
	float step3dX = (b.pos3d.x*b.zinv  - a.pos3d.x*a.zinv) / float(max(N-1,1));
	float step3dY = (b.pos3d.y*b.zinv  - a.pos3d.y*a.zinv) / float(max(N-1,1));
	float step3dZ = (b.pos3d.z*b.zinv  - a.pos3d.z*a.zinv) / float(max(N-1,1));
	float currentX = a.x, currentY = a.y, currentZ = a.zinv, 
		  current3dX = a.pos3d.x * a.zinv, current3dY = a.pos3d.y * a.zinv, current3dZ = a.pos3d.z * a.zinv;
	float minX = min(a.x,b.x);
	float maxX = max(a.x,b.x);
	float minY = min(a.y,b.y);
	float maxY = max(a.y,b.y);

	for(int i=0; i < N; ++i)
	{
		if (currentX < minX)
			currentX = minX;
		if (currentX > maxX)
			currentX = maxX;
		if (currentY < minY)
			currentY = minY;
		if (currentY > maxY)
			currentY = maxY;

		result[i].x = (int) currentX;
		result[i].y = (int) currentY;
		result[i].zinv = currentZ;
		result[i].pos3d = vec3(current3dX/currentZ, current3dY/currentZ, current3dZ/currentZ);
		currentX += stepX;
		currentY += stepY;
		currentZ += stepZ;
		current3dX += step3dX;
		current3dY += step3dY;
		current3dZ += step3dZ;
	}
}

void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels,
 vector<Pixel>& rightPixels)
{
	// 1. Find max and min y-value of the polygon
	//    and compute the number of rows it occupies.
	vector<int> y_s;
	for (int i = 0; i < vertexPixels.size(); ++i)
		y_s.push_back(vertexPixels[i].y);
	
	int minY = *min_element(y_s.begin(), y_s.end());
	int maxY = *max_element(y_s.begin(), y_s.end());

	int ROWS = maxY - minY + 1;

	// 2. Resize leftPixels and rightPixels
	//    so that they have an element for each row.
	leftPixels.resize(ROWS);
	rightPixels.resize(ROWS);

	// 3. Initialize the x-coordinates in leftPixels
	//    to some really large value and the x-coordinates
	//    in rightPixels to some really small value.

	for(int i = 0; i < ROWS; ++i) {
		leftPixels[i].x  = +numeric_limits<int>::max();
		rightPixels[i].x = -numeric_limits<int>::max();
	}

	// 4. Loop through all edges of the polygon and use
	//    linear interpolation to find the x-coordinate for
	//    each row it occupies.
	int numEdge = vertexPixels.size()*(vertexPixels.size() - 1)/2;
	vector< vector<Pixel> > edges(numEdge);
	int counter = 0;

	for(int i = 0; i < vertexPixels.size() - 1; ++i) {
		for(int j = i + 1; j < vertexPixels.size(); ++j) {
			edges[counter].resize((max(abs(vertexPixels[i].x-vertexPixels[j].x),abs(vertexPixels[i].y-vertexPixels[j].y)))+1);
			Interpolate(vertexPixels[i], vertexPixels[j], edges[counter]);
			counter++;
		}
	}

	// Update the corresponding
	// values in rightPixels and leftPixels.
	for(int i = 0; i < edges.size(); ++i) {
		for(int j = 0; j < edges[i].size(); ++j) {
			int index = edges[i][j].y - minY;
			if(leftPixels[index].x > edges[i][j].x) {
				leftPixels[index].x = edges[i][j].x;
				leftPixels[index].zinv = edges[i][j].zinv;
				leftPixels[index].pos3d = edges[i][j].pos3d;
			}

			if(rightPixels[index].x < edges[i][j].x) {
				rightPixels[index].x = edges[i][j].x;
				rightPixels[index].zinv = edges[i][j].zinv;
				rightPixels[index].pos3d = edges[i][j].pos3d;
			}
		}
	}

	// compute y-values of each point
	for(int i = 0; i < ROWS; ++i)
		leftPixels[i].y = rightPixels[i].y = i + minY;

}

// PutPixelSDL for each pixel between the start and end for each row
void DrawPolygonRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels) {
	for(int i = 0; i < leftPixels.size(); ++i) {
		// if (leftPixels[i].x == numeric_limits<int>::max() || rightPixels[i].x == -numeric_limits<int>::max())
		// 	continue;

		// multiply the whole size by 8 or so for some anti-aliasing
		vector<Pixel> row((max(abs(leftPixels[i].x-rightPixels[i].x),abs(leftPixels[i].y-rightPixels[i].y)))+1);
		Interpolate(leftPixels[i], rightPixels[i], row);

		for(const auto& point : row)
			PixelShader(point);
	}
}

// Project the vertices, compute the polygon rows, and draw them
void DrawPolygon(const vector<Vertex>& vertices) {
	int V = vertices.size();
	vector<Pixel> vertexPixels(V);
	for(int i = 0; i < V; ++i)
		VertexShader(vertices[i], vertexPixels[i]);

	vector<Pixel> leftPixels;
	vector<Pixel> rightPixels;
	ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
	DrawPolygonRows(leftPixels, rightPixels);
}


void updateCameraAngle(float angle) 
{
    yaw += angle;
    //update rotation matrix with angle
    R = mat3(vec3( cos(angle), 0, sin(angle)),
             vec3(      0,     1,   0    ),
             vec3(-sin(angle), 0, cos(angle)));
    //update camera position with rotation matrix
    cameraPos = R * cameraPos;
}