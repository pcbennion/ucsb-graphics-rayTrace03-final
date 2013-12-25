#include "Shapes.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////////////
// Read in a raw PPM file of the "P6" style.
//
// Input: "filename" is the name of the file you want to read in
// Output: "pixels" will point to an array of pixel values
//         "width" will be the width of the image
//         "height" will be the height of the image
//
// The PPM file format is:
//
//   P6
//   <image width> <image height>
//   255
//   <raw, 8-bit binary stream of RGB values>
//
// Open one in a text editor to see for yourself.
//
//////////////////////////////////////////////////////////////////////////////////
void readPPM(const char* filename, unsigned char*& pixels, int& width, int& height)
{
  // try to open the file
  FILE* file;
  file = fopen(filename, "rb");
  if (file == NULL)
  {
    cout << " Couldn't open file " << filename << "! " << endl;
    exit(1);
  }

  // read in the image dimensions
  fscanf(file, "P6\n%d %d\n255\n", &width, &height);
  int totalPixels = width * height;

  // allocate three times as many pixels since there are R,G, and B channels
  pixels = new unsigned char[3 * totalPixels];
  fread(pixels, 1, 3 * totalPixels, file);
  fclose(file);

  // output some success information
  cout << " Successfully read in " << filename << " with dimensions: "
       << width << " " << height << endl;
}

//==========================
// SPHERE CLASS
//==========================
// constructor
Sphere::Sphere(Vector3 center, float radius, Vector3 color, int material)
	{this->center = center; this->r = radius; this->color = color; this->material = material;}

// Ray Intersection
//		solves quadratic equation derived from placing ray components in the sphere equation
//		returns -1 if no intersect
//		returns # of ray lengths to closest intersect otherwise
float Sphere::findIntersect(Vector3 rayDir, Vector3 raySrc, int param)
{
	// Find vector representing center of sphere-->ray origin
	Vector3 omc = raySrc - this->center;

	// Compute three parts of quadratic equation
	float a = rayDir.dot(rayDir);
	float b = 2 * rayDir.dot(omc);
	float c = omc.dot(omc) - r*r;

	// Quadratic formula!
	//	calculate discriminant. If <0, not intersect possible
	float disc = b*b - 4*a*c;
	if(disc<0) return -1;

	// find the two solutions
	float t1 = (-b + sqrt(disc))/(2*a);
	float t2 = (-b - sqrt(disc))/(2*a);

	if(param == NEAR_SIDE)
	{
		// return smaller t
		if(t2>0 && (t2<t1 || t1<0)) return t2;
		return t1;
	} else{
		// return larger t
		if(t1>=t2) return t1;
		return t2;
	}
}

// Update function
void Sphere::updateLocation(Vector3 p, Matrix3 rot)
    {this->center = p;}

// accessors
Vector3 Sphere::getNormal(Vector3 p)
	{return Vector3(p[0] - center[0], p[1] - center[1], p[2] - center[2]);}
Vector3 Sphere::getCenter()
	{return center;}
Vector3 Sphere::getColor(Vector3 p)
    {return color;}
Vector3 Sphere::getColor()
	{return color;}
int Sphere::getMaterial()
	{return material;}


//==========================
// CYLINDER CLASS
//==========================
// constructor
Cylinder::Cylinder(Vector3 center, float radius, float length, Matrix3 rotation, Vector3 color, int material)
{
  this->center = center;
  this->r = radius; 
  this->length = length; 
  this->rotation = rotation, 
  this->color = color; 
  this->material = material;
}

// Ray Intersection
//		solves quadratic equation derived from placing ray components in the cylinder equation
//		returns -1 if no intersect
//		returns # of ray lengths to closest intersect otherwise
float Cylinder::findIntersect(Vector3 rayDir, Vector3 raySrc, int param)
{
  //convert ray to cylinder's local frame of reference
  Vector3 origin(0.0f, 0.0f, 0.0f);
  Vector3 center = this->center;
  float radius = this->r, tMin;
  Matrix3 rotation = this->rotation;
  Matrix3 rotationT = rotation.transpose();

  //convert ray to cylinder coordinates
  Vector3 raySrcCyl = rotationT * (raySrc - center);
  Vector3 rayDirCyl = rotationT * rayDir;
  Vector3 centerCyl = rotationT * center;

  //Compute intersection...
  // Compute three parts of quadratic equation
  float a = (rayDirCyl[0] * rayDirCyl[0]) + (rayDirCyl[1] * rayDirCyl[1]);
  float b = 2.0f * ((raySrcCyl[0] * rayDirCyl[0]) + (raySrcCyl[1] * rayDirCyl[1]));
  float c = (raySrcCyl[0] * raySrcCyl[0]) + (raySrcCyl[1] * raySrcCyl[1]) - (radius * radius);

  // Quadratic formula!
  //	calculate discriminant. If <0, not intersect possible
  float disc = b*b - 4*a*c;
  if(disc<0) return -1;
  float length = this->length;
  if(disc == 0){//grazed it
    tMin = -b / (2 * a);
    float z = raySrcCyl[2] + rayDirCyl[2]*tMin;
    return (z >= 0.0f && z <= length)?tMin:-1;
  }
  else{//direct hit
    float t1 = (-b - sqrt(disc)) / (2 * a);
    //Vector3 pCyl1 = raySrcCyl + rayDirCyl*t1;
    float z1 = raySrcCyl[2] + rayDirCyl[2]*t1;
    float t2 = (-b + sqrt(disc)) / (2 * a);
    //Vector3 pCyl2 = raySrcCyl + rayDirCyl*t2;
    float z2 = raySrcCyl[2] + rayDirCyl[2]*t2;
    //pick the nearest face (checking bounds)
    if(t1 <= t2 && (z1 >= 0.0f && z1 <= length)){
      tMin = t1;
    }
    else if(z2 >= 0.0f && z2 <= length){
      tMin = t2;
    }
    else{
      return -1;
    }
  }
  return tMin;
}

// Update function
void Cylinder::updateLocation(Vector3 p, Matrix3 rot)
    {this->center = p; this->rotation = rot;}

// accessors
Vector3 Cylinder::getNormal(Vector3 p){
  //...recalculate normal
  Vector3 origin(0.0f, 0.0f, 0.0f);
  Matrix3 rotation = this->rotation;
  Matrix3 rotationT = rotation.transpose();
  //place our point back in cylinder coordinates
  Vector3 pCyl = rotationT * (p - this->center);
  //Vector3 pCyl =  (rotation * p) + this->center;
  Vector3 nCyl = Vector3(pCyl[0], pCyl[1], 0.0f) - origin;
  return rotation * nCyl;    
}
Vector3 Cylinder::getCenter()
	{return center;}
Vector3 Cylinder::getColor(Vector3 p)
    {return color;}
Vector3 Cylinder::getColor()
	{return color;}
int Cylinder::getMaterial()
	{return material;}

//==========================
// TRIANGLE CLASS
//==========================
// constructor
Triangle::Triangle(Vector3 a, Vector3 b, Vector3 c, Vector3 color, int material)
{
	// save coords, color, material
	vert[0] = a;
	vert[1] = b;
	vert[2] = c;
	this->color = color;
	this->material = material;
	
	// precompute normal
	Vector3 bma(b - a);
	Vector3 cma(c - a);
	normal = bma.cross(cma);
}
// constructor for textured triangles
Triangle::Triangle(Vector3 a, Vector3 b, Vector3 c, Vector3 color, float au, float av, float bu, float bv, float cu, float cv)
{
    // save coords, color, material
    vert[0] = a;
    vert[1] = b;
    vert[2] = c;
    this->color = color;
    this->material = MAT_TEXTURE;

    // get texture stuff
    unsigned char* buf;
    readPPM("input.ppm", buf, width, height);
    texture = new float[3*width*height];
    for(int i=0;i<3*width*height;i++)
    	texture[i]=((float)(buf[i]))/255;
    this->au=au; this->av=av;
    this->bu=bu; this->bv=bv;
    this->cu=cu; this->cv=cv;

    // precompute normal
    Vector3 bma(b - a);
    Vector3 cma(c - a);
    normal = bma.cross(cma);
}

// Ray Intersection
//		Uses method found at http://geomalgorithms.com/a06-_intersect-2.html
//			to solve a simplified barycentric model
//		Returns -1 if no intersect
//		Returns # ray lengths to triangle otherwise
float Triangle::findIntersect(Vector3 rayDir, Vector3 raySrc, int param)
{
	// Calculate vectors to represent two edges of triangle
	Vector3 u = vert[1] - vert[0];
	Vector3 v = vert[2] - vert[0];

	// Find distance until intersect with triangle plane
	// 		ray-plane intersection: raySrc + r*rayDir = p, where p is a point on the plane
	Vector3 w = raySrc - vert[0]; // vector of Vert 1 --> Ray Origin
	float nw = - normal.dot(w);
	float ndir = normal.dot(rayDir);
	if(nw<0){nw = -nw; ndir = -ndir;} // ensure the normal is forward-facing
	if(abs(ndir<0.00000001)) return -1; // ray is near parallel with plane; return no intersect
	float r = nw/ndir; // distance to intersect: here's our return value
	if(r<0) return -1; // no intersect; plane is behind ray

	// Calculate useful numbers for barycentric test
	w = (raySrc + rayDir * r) - vert[0]; // Key Vertex --> Intersect Point
	float uu, uv, vv, wu, wv, denom; // dot products and denominator
	uu = u.dot(u); uv = u.dot(v); vv = v.dot(v);
	wu = w.dot(u); wv = w.dot(v);
	denom = uv*uv - uu*vv;

	// Mathematically reduced barycentric test: calculate 'edge coords' along our u and v vectors
	//		w is within triangle if s,  t, and s+t are between 0 and 1
	float s, t;
	s = (uv*wv - vv*wu)/denom;
	if(s<0||s>1) return -1;
	t = (uv*wu - uu*wv)/denom;
	if(t<0||(s+t)>1) return -1;

	// If we've made it this far, intersect is within triangle! Return the calculated distance.
	return r;
}

// Update function
void Triangle::updateLocation(Vector3 p, Matrix3 rot)
{
    // UPDATING LOCATION WOULD REQUIRE A LOT OF HARD MATH.
    // THEREFORE, THIS IS A STUB.
}

Vector3 Triangle::getColor(Vector3 p)
{
    if(material != MAT_TEXTURE)
        return color;

    float alpha, beta, gamma;
    //float x = p[0]; float y = p[1]; float z = p[2];
    Vector3 a = vert[0]; Vector3 b = vert[1]; Vector3 c = vert[2];
    int x, y;

    // get barycentric coords. Same method as intersect!
    	// Calculate vectors to represent two edges of triangle
    	Vector3 w;
    	Vector3 u = vert[1] - vert[0];
    	Vector3 v = vert[2] - vert[0];

    	// Calculate useful numbers for barycentric test
    	w = (p) - vert[0]; // Key Vertex --> Intersect Point
    	float uu, uv, vv, wu, wv, denom; // dot products and denominator
    	uu = u.dot(u); uv = u.dot(v); vv = v.dot(v);
    	wu = w.dot(u); wv = w.dot(v);
    	denom = uv*uv - uu*vv;

    	// Mathematically reduced barycentric test: calculate 'edge coords' along our u and v vectors
    	//		w is within triangle if s,  t, and s+t are between 0 and 1
    	float s, t;
    	s = (uv*wv - vv*wu)/denom;
    	t = (uv*wu - uu*wv)/denom;
    	
    	beta = s; gamma = t; alpha = 1-s-t;

   // calculate texture coords
   x = (int)((au*alpha + bu*beta + cu*gamma)*width);
   y = (int)((av*alpha + bv*beta + cv*gamma)*height);

   // set output color
   Vector3 out;
   out[0] = texture[3*width*y + 3*x + 0];
   out[1] = texture[3*width*y + 3*x + 1];
   out[2] = texture[3*width*y + 3*x + 2];
   
   //writePPM("output.ppm", texture, width, height);

   return out;
}

// accessors
Vector3 Triangle::getNormal(Vector3 p)
	{return normal;}
Vector3 Triangle::getCoord(int index)
	{return vert[index];}
Vector3 Triangle::getColor()
	{return color;}
int Triangle::getMaterial()
	{return material;}

//==========================
// LIGHT CLASS
//==========================
// constructor
Light::Light(Vector3 point, Vector3 color)
{this->point = point; this->color = color; this->length = 0; this->width = 0;}
Light::Light(Vector3 point, Vector3 color, Vector3 dirL, Vector3 dirW, int L, int W, float scale)//area light source
{this->point = point; this->color = color; this->dirLength = dirL; this->dirWidth = dirW; 
  this->length = L; this->width = W; this->scale = scale;}

// accessors
int Light::getLength(){
  return this->length;
}
int Light::getWidth(){
  return this->width;
}
Vector3 Light::getPoint(){//point light source
  return this->point;
}

Vector3 Light::getPoint(int i, int j){//area light source
  //choose random points from within each grid square
  float s = this->scale;
  return this->point + (this->dirLength * (float)i + this->dirWidth * (float)j) * s;
}
Vector3 Light::getColor()
{return color;}
