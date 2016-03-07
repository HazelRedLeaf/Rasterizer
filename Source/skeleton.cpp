#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::ivec2;

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */

const int SCREEN_WIDTH = 200;
const int SCREEN_HEIGHT = 200;
SDL_Surface* screen;
int t;
vector<Triangle> triangles;

// camera variables
vec3 cameraPos(0.f, 0.f, -3.001f);
float f = 300.f;
float cameraSpeed = 0.2f;
float yaw = -M_PI/18.f;
mat3 R;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw();
void updateCameraAngle(float angle);
void VertexShader(const vec3& v, ivec2& p);
void Interpolate(vec3 a, vec3 b, vector<vec3>& result);

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
    // if(keystate[SDLK_LEFT]) {
    //     cameraPos.x -= cameraSpeed;
    //     updateCameraAngle(-M_PI/18.f);
    // } 
    // if(keystate[SDLK_RIGHT]) {
    //     cameraPos.x += cameraSpeed;
    //     updateCameraAngle(M_PI/18.f);
    // } 
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

		for(int v=0; v<3; ++v)
		{
			ivec2 projPos;
			VertexShader(vertices[v], projPos);
			vec3 color(1,1,1);
			PutPixelSDL(screen, projPos.x, projPos.y, color);
		}	
	}

	if(SDL_MUSTLOCK(screen))
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}

void VertexShader(const vec3& v, ivec2& p)
{
	// get position of point in the camera's coordinate system
	float X = v.x - cameraPos.x;
	float Y = v.y - cameraPos.y;
	float Z = v.z - cameraPos.z;

	// project (X,Y,Z) to (x,y,f)
	p.x = f*X/Z + SCREEN_WIDTH/2;
	p.y = f*Y/Z + SCREEN_HEIGHT/2;
}

void Interpolate(vec3 a, vec3 b, vector<vec3>& result)
{
	for( unsigned int i=0; i < result.size(); i++)
	{
		result[i].x = a.x + i * (b.x - a.x) / (result.size() - 1);
		result[i].y = a.y + i * (b.y - a.y) / (result.size() - 1);
		result[i].z = a.z + i * (b.z - a.z) / (result.size() - 1);
	}
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