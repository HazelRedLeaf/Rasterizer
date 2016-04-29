#include <iostream>
#include <algorithm>
#include <glm/glm.hpp>
#include <X11/Xlib.h> 
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
vec3 cameraPos(0.f,0.f,-2.f);
float f = 300.f;
float cameraSpeed = 0.00005f;
float yaw = 0; // Yaw angle controlling camera rotation around y-axis
mat3 R;
float focus = 0.01f;

// light variables
float lightSpeed = 0.2f;
vec3 lightPos(0,-0.5,-0.7);
vec3 lightPower = 11.1f*vec3(1, 1, 1);
vec3 indirectLightPowerPerArea = 0.5f*vec3(1, 1, 1);

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw();
void DrawDepthField();
void updateCameraAngle(float angle);
void VertexShader(const Vertex& v, Pixel& p, vec3 currentNormal, vec3 currentReflectance);
void PixelShader(const Pixel& p, vec3 currentNormal, vec3 currentReflectance);
void Interpolate(Pixel a, Pixel b, vector<Pixel>& result);
void ComputePolygonRows(const vector<Pixel>& vertexPixels, 
	vector<Pixel>& leftPixels,
	vector<Pixel>& rightPixels,
	vec3 currentNormal,
	vec3 currentReflectance);
void DrawPolygonRows(const vector<Pixel>& leftPixels,
					 const vector<Pixel>& rightPixels,
					 vec3 currentNormal, vec3 currentReflectance);
void DrawPolygon(const vector<Vertex>& vertices, vec3 currentNormal, vec3 currentReflectance);
void LoadOBJModel(vector<Triangle>& triangles);

void LoadOBJModel(vector<Triangle>& triangles,string objPath,string matPath) {
    vector<vec3> tmpVs;
    vector<vec3> tmpNorms;
    vec3 mat;

    string filename = "Resources/"+objPath;
    FILE * file = fopen(filename.c_str(),"r");
    
    if(file == NULL) {
        cout << "Unable to open object file.";
        return;
    }

    while(1) {
        char type[128];
        int res = fscanf(file,"%s",type);
        if(res == EOF)
            break;
        // parse vertices
        if(strcmp(type,"v") == 0) {
            vec3 v;
            fscanf(file,"%f %f %f\n",&v.x,&v.y,&v.z);
            tmpVs.push_back(v);
        }
        // parse normals
        else if(strcmp(type,"vn")==0) {
            vec3 norm;
            fscanf(file,"%f %f %f\n",&norm.x,&norm.y,&norm.z);
            tmpNorms.push_back(norm);
        }
        // material / colour
        else if(strcmp(type,"usemtl")==0){ 
            char name[128];
            fscanf(file,"%s\n",name);
            filename = "Resources/"+matPath;
            FILE * matFile = fopen(filename.c_str(),"r");
            if(matFile == NULL) {
                cout << "Unable to open material file.";
                return;
            }
            // scan mat file
            while(1) {
                char t[128];
                res = fscanf(matFile,"%s",t);
                if(res == EOF)
                    break;
                if(strcmp(t,"newmtl")==0) {
                    fscanf(matFile,"%s\n",t);
                    // found mat
                    if(strcmp(t,name)==0) {
                        while(1) {
                            char col[128];
                            res = fscanf(matFile,"%s",col);
                            if(res == EOF)
                                break;
                            if(strcmp(col,"Kd")==0) {
                                // set mat properties
                                fscanf(matFile,"%f %f %f\n",&mat.x,&mat.y,&mat.z);
                                break;
                            }
                            
                        }
                    }
                }
            }
        }
        // parse triangle faces
        else if(strcmp(type,"f")==0) {
            string v1,v2,v3;
            unsigned int vIndex[3],normIndex[3];
            vec3 vs[3];
            vec3 norms[3];
            int matches = fscanf(file," %d//%d %d//%d %d//%d\n", &vIndex[0],&normIndex[0],&vIndex[1],&normIndex[1],&vIndex[2],&normIndex[2]);
            if(matches != 6) 
                cout << "Failed to parse a face from blender model. ";
            for(unsigned int i = 0; i < 3; i++) {
                vs[i] = tmpVs[vIndex[i]-1];
                norms[i] = tmpNorms[normIndex[i]-1];
            }
            Triangle t(vs[0],vs[1],vs[2],mat);
            // only use one norm for face
            t.normal = norms[0];
            
            triangles.push_back(t);
        }
        
    }
}

int main(int argc, char* argv[]) {
    //is necessary for multithreaded access
    XInitThreads();

	screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT);
	t = SDL_GetTicks();	// Set start value for timer.

	//Load scene triangles
    //LoadTestModel(triangles);
    LoadOBJModel(triangles, "Townhouse.obj", "Townhouse.mtl");

    //initialize camera angle with default yaw
    updateCameraAngle(yaw);

	while(NoQuitMessageSDL()) {
		Update();
		Draw();
	}

	SDL_SaveBMP( screen, "screenshot.bmp" );
	return 0;
}

void Update() {
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2-t);
	t = t2;
	cout << "Render time: " << dt << " ms." << endl;

	Uint8* keystate = SDL_GetKeyState(0);
    //Move camera
    if(keystate[SDLK_UP])
        cameraPos.z += cameraSpeed*dt;

    if(keystate[SDLK_DOWN])
        cameraPos.z -= cameraSpeed*dt;
    if(keystate[SDLK_LEFT])
        updateCameraAngle(1/dt * -M_PI/18.f);
    if(keystate[SDLK_RIGHT]) 
        updateCameraAngle(1/dt * M_PI/18.f);
    
    //Move light source
    if(keystate[SDLK_w])
        lightPos.z += lightSpeed;
    if(keystate[SDLK_s])
        lightPos.z -= lightSpeed;
    if(keystate[SDLK_d])
        lightPos.x -= lightSpeed;
    if(keystate[SDLK_a])
        lightPos.x += lightSpeed;
    if(keystate[SDLK_q])
        lightPos.y += lightSpeed;
    if(keystate[SDLK_e])
        lightPos.y -= lightSpeed;
}


void updateCameraAngle(float angle) {
    yaw += angle;
    //update rotation matrix with angle
    R = mat3(vec3( cos(angle), 0, sin(angle)),
             vec3(      0,     1,   0    ),
             vec3(-sin(angle), 0, cos(angle)));
    //update camera position with rotation matrix
    cameraPos = R * cameraPos;
}

void Draw() {
	SDL_FillRect(screen, 0, 0);

	if(SDL_MUSTLOCK(screen))
		SDL_LockSurface(screen);

	for(int y=0; y<SCREEN_HEIGHT; ++y)
		for(int x=0; x<SCREEN_WIDTH; ++x)
			depthBuffer[y][x] = 0;

	#pragma omp parallel for 
	for(size_t i = 0; i < triangles.size(); ++i) {
		vector<Vertex> vertices(3);
		vertices[0].position = triangles[i].v0;
		vertices[1].position = triangles[i].v1;
		vertices[2].position = triangles[i].v2;
		vec3 currentNormal = triangles[i].normal;
		vec3 currentReflectance = triangles[i].color;

		DrawPolygon(vertices, currentNormal, currentReflectance);
	}
    DrawDepthField();
	if(SDL_MUSTLOCK(screen))
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}

void DrawDepthField() {
    #pragma omp parallel for
    for(int y = 0; y < SCREEN_HEIGHT; ++y) {
        for(int x = 0; x < SCREEN_WIDTH; ++x) {
            if (depthBuffer[y][x] > focus - 0.001f && depthBuffer[y][x] < focus + 0.001f)
                continue;
            vec3 averageColor;
            averageColor =      GetPixelSDL(screen, x-2, y-2) +
                            4.f*GetPixelSDL(screen, x-1, y-2) +
                            7.f*GetPixelSDL(screen, x  , y-2) +
                            4.f*GetPixelSDL(screen, x+1, y-2) +
                                GetPixelSDL(screen, x+2, y-2) +

                            4.f*GetPixelSDL(screen, x-2, y-1) +
                           16.f*GetPixelSDL(screen, x-1, y-1) +
                           26.f*GetPixelSDL(screen, x  , y-1) +
                           16.f*GetPixelSDL(screen, x+1, y-1) +
                            4.f*GetPixelSDL(screen, x+2, y-1) +

                            7.f*GetPixelSDL(screen, x-2, y  ) +
                           26.f*GetPixelSDL(screen, x-1, y  ) +
                           41.f*GetPixelSDL(screen, x  , y  ) +
                           26.f*GetPixelSDL(screen, x+1, y  ) +
                            7.f*GetPixelSDL(screen, x+2, y  ) +

                            4.f*GetPixelSDL(screen, x-2, y+1) +
                           16.f*GetPixelSDL(screen, x-1, y+1) +
                           26.f*GetPixelSDL(screen, x  , y+1) +
                           16.f*GetPixelSDL(screen, x+1, y+1) +
                            4.f*GetPixelSDL(screen, x+2, y+1) +

                                GetPixelSDL(screen, x-2, y+2) +
                            4.f*GetPixelSDL(screen, x-1, y+2) +
                            7.f*GetPixelSDL(screen, x  , y+2) +
                            4.f*GetPixelSDL(screen, x+1, y+2) +
                                GetPixelSDL(screen, x+2, y+2);


            averageColor.x = averageColor.x/273;
            averageColor.y = averageColor.y/273;
            averageColor.z = averageColor.z/273;

            PutPixelSDL( screen, x, y, averageColor);
        }
    }
}

void VertexShader(const Vertex& v, Pixel& p, vec3 currentNormal, vec3 currentReflectance) {
	// get position of point in the camera's coordinate system
	float X = v.position.x - cameraPos.x;
	float Y = v.position.y - cameraPos.y;
	float Z = v.position.z - cameraPos.z;
	//vec3 P (X, Y, Z);
	vec3 Pp = vec3(X,Y,Z);

	Pp = Pp * R;

	X = Pp.x;
	Y = Pp.y;
	Z = Pp.z;

	if (Z == 0)
		return;

	// store 1/Z for the depth buffer
	p.zinv = 1/Z;
	// project (X,Y,Z) to (x,y,f)
	p.x = int(f*X/Z) + SCREEN_WIDTH/2;
	p.y = int(f*Y/Z) + SCREEN_HEIGHT/2;	
	// store the 3D position of the vertex to the corresponding
	// variable in Pixel
	p.pos3d = v.position;
}

void PixelShader(const Pixel& p, vec3 currentNormal, vec3 currentReflectance) {
	int x = p.x;
	int y = p.y;

	// the current pixel is closer to the camera than what was drawn
	// at this position before (or infinity - nothing was drawn)
	if(p.zinv > depthBuffer[y][x]) {
		//cout << "x = " << x << "\n";
		depthBuffer[y][x] = p.zinv;
		//cout << "y = " << y << "\n";

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
		//cout << "zinv = " << depthBuffer[y][x] << "\n";

		// R = rho * (D + N)
		vec3 Ref = currentReflectance * (D + indirectLightPowerPerArea);
		
		PutPixelSDL(screen, x, y, Ref);
		//cout << "After PutPixelSDL\n";
	}
}

void Interpolate_Old(Pixel a, Pixel b, vector<Pixel>& result) {
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

	for(int i=0; i < N; ++i) {
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


void Bresenham_Interpolate(Pixel a, Pixel b, vector<Pixel>& result) {
    int x1,y1,x2,y2;
    //bresenham code
    x1 = a.x;
    y1 = a.y;
    x2 = b.x;
    y2 = b.y;

    int dx = x2 - x1;
    result.resize(dx);
    int dy = y2 - y1;
    int dxdy2 = 2*dy - 2*dx;
    int d = 2 * dy - dx;
    vec3 pos3d = (b.pos3d - a.pos3d)/float(dx);

    float dzinv = (b.zinv - a.zinv)/float(dx);

    for(int i = 0; i < dx; i++) {
        x1 ++;
        if(d < 0) {
            d += 2*dy;
        }
        else {
            d += dxdy2;
            y1++;
        }
        if(x1 < SCREEN_WIDTH && x1 >= 0) {
            result[i].x = x1;
            result[i].y = y1;
            result[i].zinv = a.zinv + dzinv*float(i);
            result[i].pos3d = a.pos3d+pos3d*float(i);
        }

    }
}

void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels,
 vector<Pixel>& rightPixels, vec3 currentNormal, vec3 currentReflectance) {
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
			Interpolate_Old(vertexPixels[i], vertexPixels[j], edges[counter]);

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
void DrawPolygonRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels, vec3 currentNormal, vec3 currentReflectance) {
	for(int i = 0; i < leftPixels.size(); ++i) {
		// if (leftPixels[i].x == numeric_limits<int>::max() || rightPixels[i].x == -numeric_limits<int>::max())
		// 	continue;

		//vector<Pixel> row((max(abs(leftPixels[i].x-rightPixels[i].x),abs(leftPixels[i].y-rightPixels[i].y)))+1);
		//Interpolate_Old(leftPixels[i], rightPixels[i], row);

        vector<Pixel> row;
        Bresenham_Interpolate(leftPixels[i], rightPixels[i], row);


		//cout << "Before pixel shader\n";
		for(int j = 0; j < row.size(); j++)
			if(row[j].x >= 0 && row[j].x < SCREEN_WIDTH && row[j].y >= 0 && row[j].y < SCREEN_HEIGHT)
				PixelShader(row[j], currentNormal, currentReflectance);

	}
}

// Project the vertices, compute the polygon rows, and draw them
void DrawPolygon(const vector<Vertex>& vertices, vec3 currentNormal, vec3 currentReflectance) {
	int V = vertices.size();
	vector<Pixel> vertexPixels(V);
	for(int i = 0; i < V; ++i)
		VertexShader(vertices[i], vertexPixels[i], currentNormal, currentReflectance);

	vector<Pixel> leftPixels;
	vector<Pixel> rightPixels;
	ComputePolygonRows(vertexPixels, leftPixels, rightPixels, currentNormal, currentReflectance);
	DrawPolygonRows(leftPixels, rightPixels, currentNormal, currentReflectance);
}
