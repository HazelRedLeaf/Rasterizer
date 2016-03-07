#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::vec3;
using glm::mat3;

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */
const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
int t; // time
vector<vec3> stars(1000); 		// store location of all stars

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */
void Update();
void Draw();

int main(int argc, char* argv[])
{
	screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT);
	t = SDL_GetTicks();	// Set start value for timer (to current time)

	// create initial random positions for the stars
	for(int i = 0; i < 1000; i++)
	{
		// For x and y:
		//				float(rand()) / float(RAND_MAX) * (MAX - MIN) + MIN
		//				MAX = 1, MIN = -1
		stars[i].x = float(rand()) / float(RAND_MAX) * 2 - 1;
		stars[i].y = float(rand()) / float(RAND_MAX) * 2 - 1;
		stars[i].z = float(rand()) / float(RAND_MAX);
	}

		for(int i = 0; i < 1000; i++)
	{
		if (stars[i].x < -1 || stars[i].x > 1 ||
			stars[i].y < -1 || stars[i].y > 1 ||
			stars[i].z < 0 || stars[i].z > 1)
			cout << "( "
				 << stars[i].x << ", "
				 << stars[i].y << ", "
				 << stars[i].z << ")" << endl;
	}

	while(NoQuitMessageSDL())
	{
		Update();
		Draw();
	}

	SDL_SaveBMP(screen, "screenshot.bmp");

	return 0;
}

void Update()
{
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2-t);
	t = t2;
	//cout << "Render time: " << dt << " ms." << endl;

    // simulate motion
    for(size_t i = 0; i < stars.size(); ++i)
    {
        // update stars Z coordinate
        stars[i].z = stars[i].z - 0.005*dt;

        // a star is out of the screen
        if(stars[i].z <= 0)
            stars[i].z += 1;
        if(stars[i].z > 1)
            stars[i].z -= 1;
    }
}

void Draw()
{
	// initialise a black screen
	SDL_FillRect(screen, 0, 0);

	if(SDL_MUSTLOCK(screen))
		SDL_LockSurface(screen);

    // initial positioning of the stars
	for(size_t i = 0; i < stars.size(); ++i)
	{
		float f = SCREEN_HEIGHT / 2;	// focal length
		
		// 2D projection (u, v) of 3D coordinates (x, y, z)
		float u = f * stars[i].x / stars[i].z + SCREEN_WIDTH/ 2;
		float v = f * stars[i].y / stars[i].z + SCREEN_HEIGHT / 2;

		// position stars
		vec3 color(1.0, 1.0, 1.0);
		vec3 colorDimmed = 0.2f * vec3(1, 1, 1) / stars[i].z*stars[i].z;
		PutPixelSDL(screen, u, v, color);
	}

	if(SDL_MUSTLOCK(screen))
		SDL_UnlockSurface(screen);

	SDL_UpdateRect(screen, 0, 0, 0, 0);
}
