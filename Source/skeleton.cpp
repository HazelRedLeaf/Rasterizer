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
/* GLOBAL VARIABLES                                                            */

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
int t;
vector<Triangle> triangles;
vec3 currentColor;

// camera variables
vec3 cameraPos(0.f, 0.f, -3.001f);
float f = 500.f;
float cameraSpeed = 0.002f;
float yaw = 0; // Yaw angle controlling camera rotation around y-axis
mat3 R(vec3( 0, 0, 1),
       vec3( 0, 1, 0),
       vec3(-1, 0, 0));

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw();
void updateCameraAngle(float angle);
void VertexShader(const vec3& v, ivec2& p);
void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result);
void DrawLineSDL(SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color);
void DrawPolygonEdges(const vector<vec3>& vertices);
void ComputePolygonRows(const vector<ivec2>& vertexPixels, 
	vector<ivec2>& leftPixels,
	vector<ivec2>& rightPixels);
void DrawPolygonRows(const vector<ivec2>& leftPixels,
					 const vector<ivec2>& rightPixels);
void DrawPolygon(const vector<vec3>& vertices);

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
}

void Draw()
{
	SDL_FillRect(screen, 0, 0);

	if(SDL_MUSTLOCK(screen))
		SDL_LockSurface(screen);

	for(size_t i = 0; i < triangles.size(); ++i)
	{
		currentColor = triangles[i].color;

		vector<vec3> vertices(3);
		vertices[0] = triangles[i].v0;
		vertices[1] = triangles[i].v1;
		vertices[2] = triangles[i].v2;

		DrawPolygon(vertices);
	}

	if(SDL_MUSTLOCK(screen))
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}

void VertexShader(const vec3& v, ivec2& p)
{
	// column vectors from rotation matrix
	vec3 R1 (R[0][0], R[1][0], R[2][0]);
	vec3 R2 (R[0][1], R[1][1], R[2][1]);
	vec3 R3 (R[0][2], R[1][2], R[2][2]);

	// get position of point in the camera's coordinate system
	float X = v.x - cameraPos.x;
	float Y = v.y - cameraPos.y;
	float Z = v.z - cameraPos.z;
	//vec3 P (X, Y, Z);
	X = X * R1.x + X * R1.z;
	Z = Z * R3.x + Z * R3.z;

	// project (X,Y,Z) to (x,y,f)
	p.x = f*X/Z + SCREEN_WIDTH/2;
	p.y = f*Y/Z + SCREEN_HEIGHT/2;
}

void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result)
{
	int N = result.size();
	vec2 step = vec2(b-a) / float(max(N-1,1));
	vec2 current(a);
	float minX = min(a.x,b.x);
	float maxX = max(a.x,b.x);
	float minY = min(a.y,b.y);
	float maxY = max(a.y,b.y);

	for(int i=0; i < N; ++i)
	{
		if (current.x < minX)
			current.x = minX;
		if (current.x > maxX)
			current.x = maxX;
		if (current.y < minY)
			current.y = minY;
		if (current.y > maxY)
			current.y = maxY;

		result[i] = current;
		current += step;
	}
}

void DrawLineSDL(SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color)
{
	ivec2 delta = glm::abs(a-b);
	int pixels = glm::max(delta.x, delta.y) + 1;
	vector<ivec2> line(pixels);
	Interpolate(a, b, line);

	for (int i = 0; i < pixels; ++i)
		PutPixelSDL(screen, line[i].x, line[i].y, color);
}

void DrawPolygonEdges(const vector<vec3>& vertices)
{
	int V = vertices.size();
	// Transform each vertex from 3D world position to 2D image position:
	vector<ivec2> projectedVertices(V);
	for(int i = 0; i < V; ++i)
	{
		VertexShader(vertices[i], projectedVertices[i]);
	}

	// Loop over all vertices and draw the edge from it to the next vertex:
	for(int i = 0; i < V; ++i)
	{
		int j = (i + 1) % V; // The next vertex
		vec3 color( 1, 1, 1 );
		DrawLineSDL(screen, projectedVertices[i], projectedVertices[j], color);
	}
}

void ComputePolygonRows(const vector<ivec2>& vertexPixels, vector<ivec2>& leftPixels,
 vector<ivec2>& rightPixels)
{
	// 1. Find max and min y-value of the polygon
	//    and compute the number of rows it occupies.
	vector<int> x_s,y_s;
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
	//    each row it occupies. Update the corresponding
	//    values in rightPixels and leftPixels.
	vector< vector<ivec2> > edges(vertexPixels.size()*(vertexPixels.size() - 1)/2);
	int counter = 0;

	for(int i = 0; i < vertexPixels.size() - 1; ++i) {
		for(int j = i + 1; j < vertexPixels.size(); ++j) {
			edges[counter].resize((max(abs(vertexPixels[i].x-vertexPixels[j].x),abs(vertexPixels[i].y-vertexPixels[j].y)))+1);
			Interpolate(vertexPixels[i], vertexPixels[j], edges[counter]);
			counter++;
		}
	}

	for(int i = 0; i < edges.size(); ++i) {
		for(int j = 0; j < edges[i].size(); ++j) {
			int index = edges[i][j].y - minY;
			leftPixels[index].x = min(leftPixels[index].x , edges[i][j].x);
			rightPixels[index].x = max(rightPixels[index].x , edges[i][j].x);
		}
	}

	for(int i = 0; i < ROWS; ++i)
		leftPixels[i].y = rightPixels[i].y = i + minY;
}

// PutPixelSDL for each pixel between the start and end for each row
void DrawPolygonRows(const vector<ivec2>& leftPixels, const vector<ivec2>& rightPixels) {
	for(int i = 0; i < leftPixels.size(); ++i) {
		vector<ivec2> row((max(abs(leftPixels[i].x-rightPixels[i].x),abs(leftPixels[i].y-rightPixels[i].y)))+1);
		Interpolate(leftPixels[i], rightPixels[i], row);

		for(const ivec2& point : row)
			PutPixelSDL(screen, point.x, point.y, currentColor);
	}
}

void DrawPolygon(const vector<vec3>& vertices) {
	int V = vertices.size();
	vector<ivec2> vertexPixels(V);
	for(int i = 0; i < V; ++i)
		VertexShader(vertices[i], vertexPixels[i]);

	vector<ivec2> leftPixels;
	vector<ivec2> rightPixels;
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