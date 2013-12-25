//////////////////////////////////////////////////////////////////////////////////
// This is a front end for a set of viewer clases for the Carnegie Mellon
// Motion Capture Database: 
//    
//    http://mocap.cs.cmu.edu/
//
// The original viewer code was downloaded from:
//
//   http://graphics.cs.cmu.edu/software/mocapPlayer.zip
//
// where it is credited to James McCann (Adobe), Jernej Barbic (USC),
// and Yili Zhao (USC). There are also comments in it that suggest
// and Alla Safonova (UPenn) and Kiran Bhat (ILM) also had a hand in writing it.
//
// I have made minor modifications throughout so that you can get direct access 
// to the bone transforms in the "display()" function below.
//
// Professor Kim
// kim@mat.ucsb.edu
// 11/15/13
//
//////////////////////////////////////////////////////////////////////////////////
//osx
/*#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>*/

//linux
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>






#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <ctime>
//#include <Windows.h> // for some reason, this is where time(NULL) is O_o

#include "Vector3.h"
#include "Matrix3f.h"
#include "Shapes.h"

// Dr.Kim's struff
#include "MATRIX4.h"
#include "MATRIX3.h"
#include "skeleton.h"
#include "displaySkeleton.h"
#include "motion.h"

using namespace std;

GLfloat cyan[] = { 0.0, 1.0, 1.0, 1.0 };
GLfloat magenta[] = { 1.0, 0.0, 1.0, 1.0 };
GLfloat white[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat red[] = { 1, 0, 0, 1.0 };
GLfloat blue[] = { 0.0, 0.0, 1.0, 1.0 };
GLfloat green[] = { 0.0, 1.0, 0.0, 1.0 };

// viewport variables
unsigned char* pixels;
int width, height;

// camera info
Vector3 eye; Vector3 focus; Vector3 up;
float fovy;
float nearz, nearx, neary;

// scene objects
Shape* objects[500];
Cylinder* cylinders[100];
Sphere* spheres[100];
Triangle* triangles[100];
Light* lights[100];
int numObjects;
int numSpheres;
int numCylinders;
int numTriangles;
int numLights;

// Stick-man classes
DisplaySkeleton displayer;		
Skeleton* skeleton;
Motion* motion;

// bone animation stuff
int boneIndex;
int objsPerBone;

int windowWidth = 640;
int windowHeight = 480;

// cameraAngle of rotation for the camera direction
float cameraAngle = 1.57;
// actual vector representing the camera's direction
float lx = 1.0f, lz = 0.0f;
// XZ position of the camera
float xPos = -6.0f, zPos = 1.0f;

GLUquadricObj *obj;

// current motion capture frame that is being displayed
int currentFrameIndex = 0;
int maxFrameIndex = 1800;
bool animate = true;

#define PI 3.14159265
static float toRadians = M_PI/180.0f;

//////////////////////////////////////////////////////////////////////////////////
////Writes a new PPM file from the current frame buffer
//////////////////////////////////////////////////////////////////////////////////
void writePPM()
{
  static int frameNumber = 0;

  char buffer[256];
  sprintf(buffer, "./frames/frame.%04i.ppm", frameNumber);
  cout << "Writing file " << buffer << " ... " << flush;

   //allocate the image array
  FILE* file;
  file = fopen(buffer, "wb");
  if (file == NULL)
    cout << " File write failed! " << endl;

  fprintf(file, "P6\n%d %d\n255\n", windowWidth, windowHeight);

  int totalPixels = windowWidth * windowHeight;

  GLubyte* frameBuffer = new GLubyte[totalPixels * 3];
  glReadPixels(0, 0, windowWidth, windowHeight, GL_RGB,GL_UNSIGNED_BYTE, frameBuffer);
  unsigned char* image = new unsigned char[3 * totalPixels];
  for (int y = 0; y < windowHeight; y++)
    for (int x = 0; x < windowWidth; x++)
    {
      int index = x + y * windowWidth;
      int reverse = x + ((windowHeight - 1 - y) * windowWidth);

      image[3 * reverse] = frameBuffer[3 * index];
      image[3 * reverse + 1] = frameBuffer[3 * index + 1];
      image[3 * reverse + 2] = frameBuffer[3 * index + 2];
    }
  fwrite(image, 1, 3 * totalPixels, file);
  fclose(file);
  frameNumber++;

  delete[] image;
  cout << " done." << endl;
}


//////////////////////////////////////////////////////////////////////////////////
//// Compute the current modelview matrix
//////////////////////////////////////////////////////////////////////////////////
Matrix3 getModelviewMatrix()
{
    up /= up.length(); // ensure the up vector is normalized

    // compute x, y, and z of modelview transform
    Vector3 foc = focus; //foc[1]*=-1;
    Vector3 f = (foc)-eye;
        f/=f.length();
    Vector3 s = f.cross(up);
    Vector3 u = (s/s.length()).cross(f);

    // correct dimensions
    f = -f;

    return Matrix3(s, u, f);
}

//////////////////////////////////////////////////////////////
// Return a 3x3 rotation matrix
//////////////////////////////////////////////////////////////
Matrix3 getRotationMatrix(float angle, float X, float Y, float Z){
  Vector3 axis(X, Y, Z);
  axis/=axis.length();

  float x = axis[0], y = axis[1], z = axis[2];//normalized x, y, z
  float c = cos(toRadians * angle), s = sin(toRadians * angle);

  //create rotation matrix
  float r[9] = {
    (x*x)*(1-c)+c,     (y*x)*(1-c)+(z*s), (x*z)*(1-c)-(y*s), //col 1
    (x*y)*(1-c)-(z*s), (y*y)*(1-c)+c,     (y*z)*(1-c)+(x*s), //col 2
    (x*z)*(1-c)+(y*s), (y*z)*(1-c)-(x*s), (z*z)*(1-c)+c     //col 3
  };
  return Matrix3(r);
}

//////////////////////////////////////////////////////////////////////////////////
//// Helper functions for raytracing
//////////////////////////////////////////////////////////////////////////////////
// find closest object that intersects with ray
float findClosestIntersect(Vector3 rayDir, Vector3 raySrc, int &object)
{
    float t = 1000, current;
    for(int i=0; i<numObjects; i++)
    {
        current = objects[i]->findIntersect(rayDir, raySrc, NEAR_SIDE);
        if(current > .01 && current < t)
        {
            // get return value and object index
            t = current;
            object = i;
        }
    }
    if(t>=999) return -1;
    return t;
}
// Shading based on Blinn-Phong lighting model, as described in the book and in class
//		src is the ray source, p is the point in question, ints are the indexes of object and light
//		Returns the amount that the color vector should be incremented
Vector3 calcShading(Vector3 p, Vector3 src, Shape* object, Light* light)
{
	Vector3 out = Vector3(0, 0, 0);

	// Start by calculating a bunch of useful normalized vectors
	Vector3 n = object->getNormal(p); n /= n.length(); // normal at p
	Vector3 v = src - p; v /= v.length(); // ray P-->Eye

	//get dimensions of light
	int length = light->getLength(), width = light->getWidth(); 
	if(length * width == 0){// point light calculations
	  Vector3 l = light->getPoint() - p; l /= l.length(); // ray P-->Light Source
	  
	  // Make sure normal faces forward
	  if(v.dot(n)<0) n = -n;
	  
	  // Calculate "perfect reflection" vector and some cosines
	  float nl = n.dot(l); nl = (nl>0) ? nl : 0;
	  Vector3 r = n*nl*2 - l; r /= r.length();
	  float rv = r.dot(v);
	  
	  // If material is not diffuse, zero out diffuse component
	  if(object->getMaterial() > MAT_DIFFUSE) nl = 0;
	  
	  // If there's nothing between p and the light, add shading
	  int o;
	  if(findClosestIntersect((p-light->getPoint()), light->getPoint(), o)>.99)
	    {
	      Vector3 blend = object->getColor(p) * light->getColor();
	      // Increment output color. Assumes Kd and Ks are both 1
	      out += /* Diffuse: */	blend*nl
		+  /* Specular:*/	blend*pow(rv,10);
	    }
	  return out;
	}
	else{// distributed area light
	  int raysPer = 20;
	  int numRays = length * width * raysPer;
	  float jitter = 1.0f / 10000.0f; 
	  //iterate over the light's grid
	  for(int i = 0; i < length; i++){
	    for(int j = 0; j < width; j++){
	      int k = 0;
	      while(k < raysPer){	 
		//get random offset
		srand(time(NULL) + rand());
		float offsetX = (float)((rand() % 500)) * jitter;
		float offsetY = (float)((rand() % 500)) * jitter;
		float offsetZ = (float)((rand() % 500)) * jitter;
		
		Vector3 lp = light->getPoint(i, j) + Vector3(offsetX, offsetY, offsetZ);;
		//cout << "offsets: {" << offsetX << ", " << offsetY << ", " << offsetZ << "}" << endl; 
		//cout << "light point " << light->getPoint() << " at grid " << i << ", " << j << ": " << lp << endl;

		Vector3 l = lp - p; l /= l.length(); // ray P-->Light Source

		// Make sure normal faces forward
		if(v.dot(n)<0) n = -n;
		
		// Calculate "perfect reflection" vector and some cosines
		float nl = n.dot(l); nl = (nl>0) ? nl : 0;
		Vector3 r = n*nl*2 - l; r /= r.length();
		float rv = r.dot(v);

		// If material is not diffuse, zero out diffuse component
		if(object->getMaterial() > MAT_DIFFUSE) nl = 0;
		
		// If there's nothing between p and the light, add shading
		int o;
		
		if(findClosestIntersect((p - lp), lp, o)>.99)
		  {
		    Vector3 blend = object->getColor(p) * light->getColor();
		    // Increment output color. Assumes Kd and Ks are both 1
		    out += /* Diffuse: */	blend*nl
		      +  /* Specular:*/	blend*pow(rv,10);
		  }
		k++;
	      }
	    }
	  }
	  return out /= numRays;
	}
}
// Reflection & Refraction
Vector3 reflectRefract(Vector3 p, Vector3 rayDir, Vector3 raySrc, Shape* object, int depth)
{
	// if past max depth, color dark grey
	if(depth>9) return Vector3(0.2, 0.2, 0.2);

	Vector3 reflectColor = Vector3(0, 0, 0);
	Vector3 refractColor = Vector3(0, 0, 0);
	Vector3 nextP, t;
	Shape* nextObj;
	float dist;
	int i, o;

	// REFLECTION:
	// find reflected ray
	Vector3 n = object->getNormal(p); n /= n.length();
	if(n.dot(rayDir)<0) n = -n;
	rayDir /= rayDir.length();
	float dn = n.dot(rayDir);
	Vector3 r = rayDir-n*dn*2; r /= r.length();
	//generate jittered reflection rays
	int k = 0, raysPer = 30;
	float jitter = 1.0f / 10000.0f; 
	while(k < raysPer){
	  srand(time(NULL) + rand());
	  float offsetX = (float)((rand() % 1000) - 500) * jitter;
	  float offsetY = (float)((rand() % 1000) - 500) * jitter;
	  float offsetZ = (float)((rand() % 1000) - 500) * jitter;
	  //cout << "offsets: {" << offsetX << ", " << offsetY << ", " << offsetZ << "}" << endl; 
	  Vector3 rJit = Vector3(r[0]+offsetX, r[1]+offsetY, r[2]+offsetZ);
			  rJit/=rJit.length();
	  if(rJit.dot(n) < 0.0f){
	    // find closest object along reflected ray
	    dist = findClosestIntersect(rJit, p, o);
	    if(dist>0)
	      {
		// if there is an intersect, get that point
		nextP = p + rJit*dist;
		nextObj = objects[o];
		// if the object is diffuse, calc shading for all lights
		if(nextObj->getMaterial()<=MAT_DIFFUSE)
		  {for(i=0; i<numLights; i++) reflectColor += calcShading(nextP, p, nextObj, lights[i]);}
		// otherwise, reflect/refract again!
		else reflectColor += reflectRefract(nextP, rJit, p, nextObj, depth+1);
	      }
	    else{
	      o = 0;
	    }
	    k++;
	  }
	}
	reflectColor /= raysPer;
	// REFRACTION:
	float mixCoeff, root;
	float nr = 1 / 1.5;
	float R0 = ((nr-1)/(nr+1))*((nr-1)/(nr+1));
	rayDir = p - raySrc; rayDir /= rayDir.length();
	n = object->getNormal(p); n /= n.length();
	dn = n.dot(rayDir);

	if(object->getMaterial()==MAT_REFRACT && dn<0)
	{
		// get refracted ray
		dn = -dn;
        root = sqrt(1 - (nr*nr)*(1 - (dn*dn)));
		if(root>0)
		{
            t = rayDir*nr - n*(nr*dn + root); t /= t.length();
			// now find other side of object
			dist = object->findIntersect(t, p, FAR_SIDE);
			nextP = p + t*dist;
			refractColor += reflectRefract(nextP, t, p, object, depth+1);
		}
		// calculate mixing coefficient for fresnel effects
		mixCoeff = R0 + (1-R0)*(1-dn)*(1-dn)*(1-dn)*(1-dn)*(1-dn);
	}else if(object->getMaterial()==MAT_REFRACT && dn>0)
	{
		// get refracted ray
        nr = 1/nr; n = -n;// reverse quantities since we're exiting
        root = sqrt(1 - (nr*nr)*(1 - (dn*dn)));
        if(root>0)
        {
            t = rayDir*nr - n*(nr*dn + root); t /= t.length();
			// now shoot t out of p to find next object
			dist = findClosestIntersect(t, p, o);
			if(dist>0)
			{
				nextP = p + t*dist;
				nextObj = objects[o];
				// if the object is diffuse, calc shading for all lights
				if(nextObj->getMaterial()<=MAT_DIFFUSE)
					{for(int i=0; i<numLights; i++) refractColor += calcShading(nextP, p, nextObj, lights[i]);}
				// otherwise, reflect/refract again!
				else refractColor += reflectRefract(nextP, rayDir, p, nextObj, depth+1);
			}
		}
		// calculate mixing coefficient for fresnel effects
		mixCoeff = R0 + (1-R0)*(1-dn)*(1-dn)*(1-dn)*(1-dn)*(1-dn);
	} else mixCoeff = 1;

	// blend reflect color and refract color, plus the normal shading for each light
	Vector3 out = Vector3(0, 0, 0);
	Vector3 shading = Vector3(0, 0, 0);
	for(i=0; i<numLights; i++) shading += calcShading(p, raySrc, object, lights[i]);
	out = shading + (reflectColor*(mixCoeff) + refractColor*(1-mixCoeff));

	return out;
}

//////////////////////////////////////////////////////////////////////////////////
// Draws to the OpenGL window
//////////////////////////////////////////////////////////////////////////////////
void display()
{
    // precompute all the bone positions of the skeleton
    displayer.ComputeBonePositionsManual(DisplaySkeleton::BONES_AND_LOCAL_FRAMES);

    // retrieve all the bones of the skeleton
    vector<MATRIX4>& rotations  = displayer.rotations();
    Matrix3 rotation;
    vector<MATRIX4>& scalings   = displayer.scalings();
    Matrix3 scaling;
    vector<VEC4F>& translations = displayer.translations();
    Vector3 translation;
    vector<float>& lengths      = displayer.lengths();
    int totalBones = rotations.size();
    objsPerBone = 2;
    // if bones are uninitialized, run inits
    int i, obj;
    if(boneIndex<0)
    {
      boneIndex=numObjects;
      for(i=1;i<totalBones;i++)
      {
          rotation = rotations[i];
          scaling = scalings[i];
          translation = translations[i];
	  if(i == 23 || i == 30){//draw boxing gloves
	    spheres[numSpheres] = new Sphere(translation, 0.075, Vector3(1, 0, 0), MAT_DIFFUSE);
              objects[numObjects++] = spheres[numSpheres++];
	    cylinders[numCylinders] = new Cylinder(translation, 0.0001, 0.0001, rotation,Vector3( 1,0,0),MAT_DIFFUSE);
              objects[numObjects++] = cylinders[numCylinders++];	      
	  }
	  else{
	    spheres[numSpheres] = new Sphere(translation, 0.025, Vector3(1, 0, 0), MAT_DIFFUSE);
              objects[numObjects++] = spheres[numSpheres++];
	    cylinders[numCylinders] = new Cylinder(translation,0.025,lengths[i],rotation,Vector3( 1,0,0),MAT_DIFFUSE);
              objects[numObjects++] = cylinders[numCylinders++];
	  }
      }
      // update bone positions
    } else {
        for(i=1;i<totalBones;i++)
        {
            rotation = rotations[i];
            scaling = scalings[i];
            translation = translations[i];
            for(obj=0;obj<objsPerBone;obj++)
            {
                objects[boneIndex + (i-1)*objsPerBone + obj]->updateLocation(translation, rotation);
            }
        }
    }

    // move camera
    eye = Vector3(0,5,-1);
    focus = Vector3(0, 3, -1);
    up = Vector3(0, 1, 0);
    rotation = rotations[16]; scaling = scalings[16]; translation = translations[16];
    eye = (rotation*scaling) * eye + (translation);
    focus = (rotation*scaling) * focus + (translation);
    focus[1]=eye[1];
    //up = (rotation*scaling) * up + (translation);

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  int index;
  int object;

  Vector3 rayDir, color, p, n, tmp;
  Matrix3 modelview = getModelviewMatrix();

  float dist;

  float r = nearx, l = -nearx;
  float t = neary, b = -neary;
  float u, v;
  int x, y;
  for(x=0; x<width; x++)
  {
      for(y=0; y<height; y++)
      {
          index = 3*width*y + 3*x;

          // calculate ray direction
          u = l + (r-l)*(x+0.5)/width;
          v = b + (t-b)*(y+0.5)/height;
          rayDir[0] = u; rayDir[1] = v; rayDir[2] = -nearz;
          rayDir = modelview.transpose()*rayDir;

          // find closest object intersect point
          dist = findClosestIntersect(rayDir, eye, object);

          // if intersect exists, determine pixel color from lighting
          if(dist>=1)
          {
              // calculate coords of collision, normal vector
              p = eye + rayDir * dist;
              n = objects[object]->getNormal(p);

              // coloring computations
              color[0] = 0; color[1] = 0; color[2] = 0; // zero out pixel color
              // if there is ambient light, it goes in this line
              
	      //FOR RENDERING WITHOUT EFFECTS:
              //color = objects[object]->getColor(p);

              // if the object is diffuse, calc shading for all lights
              if(objects[object]->getMaterial()<=MAT_DIFFUSE)
                  {for(int i=0; i<numLights; i++) color += calcShading(p, eye, objects[object], lights[i]);}
              // otherwise, reflect/refract!
              else color += reflectRefract(p, rayDir, eye, objects[object], 0);

              // prevent color components from going out of bounds
              color[0] = (color[0]>1) ? 1 : color[0];
              color[1] = (color[1]>1) ? 1 : color[1];
              color[2] = (color[2]>1) ? 1 : color[2];

              // place pixel colors in pixel buffer
              pixels[index+0]=(int)(color[0]*255);
              pixels[index+1]=(int)(color[1]*255);
              pixels[index+2]=(int)(color[2]*255);

          } else {pixels[index+0]=0; pixels[index+1]=0; pixels[index+2]=0;} // if no intersect, color black
      }
  }

  // Render picture
  glDrawPixels(width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels);

  glutSwapBuffers();
  writePPM();
}

//////////////////////////////////////////////////////////////////////////////////
// Handles keyboard events
//////////////////////////////////////////////////////////////////////////////////
void keyboard(unsigned char k, int x, int y)
{
  switch (k)
  {
    // the escape key and 'q' quit the program
    case 27:
    case 'q':
     //movie.writeMovie("stickman.mov");
      exit(0);
      break;
//    case 'a':
//      animate = !animate;
//      break;
  }
  glutPostRedisplay();
}

//////////////////////////////////////////////////////////////////////////////////
// Load up a new motion captured frame
//////////////////////////////////////////////////////////////////////////////////
void SetSkeletonsToSpecifiedFrame(int frameIndex)
{
  if (frameIndex < 0)
  {
    printf("Error in SetSkeletonsToSpecifiedFrame: frameIndex %d is illegal.\n", frameIndex);
    exit(0);
  }
  if (displayer.GetSkeletonMotion(0) != NULL)
  {
    int postureID;
    if (frameIndex >= displayer.GetSkeletonMotion(0)->GetNumFrames())
      postureID = displayer.GetSkeletonMotion(0)->GetNumFrames() - 1;
    else 
      postureID = frameIndex;
    displayer.GetSkeleton(0)->setPosture(* (displayer.GetSkeletonMotion(0)->GetPosture(postureID)));
  }
}

//////////////////////////////////////////////////////////////////////////////////
// Called occasionally to see if anything's happening
//////////////////////////////////////////////////////////////////////////////////
void idle()
{
  // load up the next motion capture frame
  SetSkeletonsToSpecifiedFrame(currentFrameIndex);

  if (animate)
    currentFrameIndex+=2;

  // if we've reached the end of the motion capture sequence,
  // start over from frome 0
  int totalFrames = displayer.GetSkeletonMotion(0)->GetNumFrames();
  if (currentFrameIndex >= maxFrameIndex)
  {
    currentFrameIndex = 0;
    cout << " Hit end frame, exiting" << endl;
    exit(0);
  }

  glutPostRedisplay();
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
    // init frame dimensions and buffer
    width = windowWidth;//800;
    height = windowHeight;//600;
    pixels = new unsigned char[3*width*height];

    if(argc==3)
    {
        currentFrameIndex = atoi(argv[1]);
        maxFrameIndex = atoi(argv[2]);
    }


    // init camera
    fovy = 65;
    nearz = 1;
    neary = nearz * tan(fovy/2 * (PI/180));
    nearx = neary * ((float)width/(float)height);
    eye = Vector3(-5, 1, -5);
    focus= Vector3(0, 1, 0);
    up = Vector3(0, 1, 0);

    // init scene objects
    numObjects = 0;
    boneIndex=-1;
    objsPerBone=0;
    // spheres
    numSpheres=0;
    spheres[numSpheres] = new Sphere(Vector3(0, -1000, 0), 1000, Vector3(0.5, 0.5, 0.5), MAT_DIFFUSE);
    objects[numObjects++] = spheres[numSpheres++];
    spheres[numSpheres] = new Sphere(Vector3(-1.5, 2, 3), 0.2, Vector3(0, 0, 1), MAT_DIFFUSE);
    objects[numObjects++] = spheres[numSpheres++];
    //cylinders
    numCylinders=0;
    Matrix3 rot = getRotationMatrix(90, 1, 0, 0);
    cylinders[numCylinders]= new Cylinder(Vector3(-1.5, 0.5, 3), 0.2, 1.5, rot, Vector3(0, 0, 1), MAT_DIFFUSE);
    objects[numObjects++] = cylinders[numCylinders++];
    cylinders[numCylinders]= new Cylinder(Vector3(-1.5, 2, 3), 0.005, 2, rot, Vector3(0.8, 0.8, 0.8), MAT_DIFFUSE);
    objects[numObjects++] = cylinders[numCylinders++];
    // dumbells
    float l = 0.35;
    Vector3 pos = Vector3(-3, 0.07, 3);
    rot = getRotationMatrix(120, 0, 1, 0);
    spheres[numSpheres] = new Sphere(pos, 0.07, Vector3(0.2, 0.2, 0.2), MAT_DIFFUSE);
    objects[numObjects++] = spheres[numSpheres++];
    spheres[numSpheres] = new Sphere(pos+(rot*Vector3(0, 0, l)), 0.07, Vector3(0.2, 0.2, 0.2), MAT_DIFFUSE);
    objects[numObjects++] = spheres[numSpheres++];
    cylinders[numCylinders]= new Cylinder(pos, 0.02, l, rot, Vector3(0.2, 0.2, 0.2), MAT_DIFFUSE);
    objects[numObjects++] = cylinders[numCylinders++];
    
    pos = Vector3(-3.5, 0.07, 3.4);
    rot = getRotationMatrix(170, 0, 1, 0);
    spheres[numSpheres] = new Sphere(pos, 0.07, Vector3(0.2, 0.2, 0.2), MAT_DIFFUSE);
    objects[numObjects++] = spheres[numSpheres++];
    spheres[numSpheres] = new Sphere(pos+(rot*Vector3(0, 0, l)), 0.07, Vector3(0.2, 0.2, 0.2), MAT_DIFFUSE);
    objects[numObjects++] = spheres[numSpheres++];
    cylinders[numCylinders]= new Cylinder(pos, 0.02, l, rot, Vector3(0.2, 0.2, 0.2), MAT_DIFFUSE);
    objects[numObjects++] = cylinders[numCylinders++];

    //triangles
    numTriangles=0;
    Vector3 ah = Vector3( 2, 3,  5); Vector3 al = Vector3( 2, 0,  5);
    Vector3 bh = Vector3(-5, 3,  5); Vector3 bl = Vector3(-5, 0,  5);
    Vector3 ch = Vector3( 2, 3, -2); Vector3 cl = Vector3( 2, 0, -2);
    Vector3 dh = Vector3(-5, 3, -2); Vector3 dl = Vector3(-5, 0, -2);
    triangles[numTriangles] = new Triangle(ah, bh, al, Vector3(0, 0, 1), 0, 0, 1, 0, 0, 1);
    objects[numObjects++] = triangles[numTriangles++];
    triangles[numTriangles] = new Triangle(bh, al, bl, Vector3(0, 0, 1), 1, 0, 0, 1, 1, 1);
    objects[numObjects++] = triangles[numTriangles++];
    triangles[numTriangles] = new Triangle(bh, dh, bl, Vector3(0.5, 0.5, 0.5),MAT_REFLECT);
    objects[numObjects++] = triangles[numTriangles++];
    triangles[numTriangles] = new Triangle(dh, bl, dl, Vector3(0.5, 0.5, 0.5),MAT_REFLECT);
    objects[numObjects++] = triangles[numTriangles++];
    triangles[numTriangles] = new Triangle(ch, ah, cl, Vector3(1, 0, 1), 0, 0, 1, 0, 0, 1);
    objects[numObjects++] = triangles[numTriangles++];
    triangles[numTriangles] = new Triangle(ah, cl, al, Vector3(1, 0, 1), 1, 0, 0, 1, 1, 1);
    objects[numObjects++] = triangles[numTriangles++];
    triangles[numTriangles] = new Triangle(dh, ch, dl, Vector3(1, 1, 0), 0, 0, 1, 0, 0, 1);
    objects[numObjects++] = triangles[numTriangles++];
    triangles[numTriangles] = new Triangle(ch, dl, cl, Vector3(1, 1, 0), 1, 0, 0, 1, 1, 1);
    objects[numObjects++] = triangles[numTriangles++];
    triangles[numTriangles] = new Triangle(ah, bh, ch, Vector3(0.75, 0.75, 0.75),MAT_DIFFUSE);
    objects[numObjects++] = triangles[numTriangles++];
    triangles[numTriangles] = new Triangle(bh, ch, dh, Vector3(0.75, 0.75, 0.75),MAT_DIFFUSE);
    objects[numObjects++] = triangles[numTriangles++];
    // lights
    numLights=0;
    lights[numLights++] = new Light(Vector3(1, 2, -1), Vector3(0.5, 0.5, 0.5), 
				    Vector3(1.0f, 0.0f, 0.0f), Vector3(0.0f, 0.0f, 1.0f), 1, 1, 0.05f);
    lights[numLights++] = new Light(Vector3(-1, 2.9, -1), Vector3(0.5, 0.5, 0.5), 
				    Vector3(1.0f, 0.0f, 0.0f), Vector3(0.0f, 0.0f, 1.0f), 1, 1, 0.05f);
    
    // GLUT functions
  glutInit(&argc, argv);
  glutInitWindowSize(windowWidth,windowHeight);
  glutInitWindowPosition(100, 100);
  glutCreateWindow("CS 180, Final Project");
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
  glutDisplayFunc(display);
  glutKeyboardFunc(keyboard);
  glutIdleFunc(idle);
  
  // load up Kim's skeleton stuff
  string skeletonFilename("01.asf");
  string motionFilename("01_02.amc");
  skeleton = new Skeleton(skeletonFilename.c_str(), MOCAP_SCALE);
  skeleton->setBasePosture();
  displayer.LoadSkeleton(skeleton);
  motion = new Motion(motionFilename.c_str(), MOCAP_SCALE, skeleton);
  displayer.LoadMotion(motion);
  skeleton->setPosture(*(displayer.GetSkeletonMotion(0)->GetPosture(0)));
  // load up the next motion capture frame
  SetSkeletonsToSpecifiedFrame(currentFrameIndex);

  glutMainLoop();

  return 0;
}
