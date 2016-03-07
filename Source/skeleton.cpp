#include <iostream>
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

// camera variables
vec3 cameraPos(0.32f, 0.f, -1.8f);
float f = 200.f;
float cameraSpeed = 0.2f;
float yaw = -M_PI/18.f;
mat3 R(vec3(1, 1, 1),
       vec3(1, 1, 1),
       vec3(1, 1, 1));

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw();
void updateCameraAngle(float angle);
void VertexShader(const vec3& v, ivec2& p);
void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result);
void DrawLineSDL(SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color);

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
		vector<vec3> vertices(3);
		vertices[0] = triangles[i].v0;
		vertices[1] = triangles[i].v1;
		vertices[2] = triangles[i].v2;
		vector<ivec2> projPoss;

		vec3 color(1,1,1);

		for(int v=0; v<3; ++v)
		{
			ivec2 projPos;
			VertexShader(vertices[v], projPos);
			projPoss.push_back(projPos);
			PutPixelSDL(screen, projPos.x, projPos.y, color);
		}

		DrawLineSDL(screen, projPoss[0], projPoss[1], color);
		DrawLineSDL(screen, projPoss[0], projPoss[2], color);
		DrawLineSDL(screen, projPoss[1], projPoss[2], color);	
		projPoss.clear();
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
	// X = X * R1.x + Z * R1.z;
	// Z = X * R3.x + Z * R3.z;

	// project (X,Y,Z) to (x,y,f)
	p.x = f*X/Z + SCREEN_WIDTH/2;
	p.y = f*Y/Z + SCREEN_HEIGHT/2;
}

void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result)
{
	int N = result.size();
	vec2 step = vec2(b-a) / float(max(N-1,1));
	vec2 current(a);

	for(int i=0; i < N; ++i)
	{
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

	for (size_t i = 0; i < pixels; i++)
		PutPixelSDL(screen, line[i].x, line[i].y, color);
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