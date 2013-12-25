/*
Revision 1 - Steve Lin, Jan. 14, 2002
Revision 2 - Alla and Kiran, Jan 18, 2002
Revision 3 - Jernej Barbic and Yili Zhao, Feb, 2012
*/
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include "types.h"


//osx
/*#include <OpenGL/gl.h>
  #include <GLUT/glut.h>*/

//linux
#include <GL/gl.h>
#include <GL/glut.h>


#include "skeleton.h"
#include "motion.h"
#include "displaySkeleton.h"

float DisplaySkeleton::jointColors[NUMBER_JOINT_COLORS][3] =
{
  {0.0f, 1.0f, 0.0f},  // GREEN
  {1.0f, 0.0f, 0.0f},  // RED
  {0.0f, 0.0f, 1.0f}   // BLUE
};

DisplaySkeleton::DisplaySkeleton(void)
{
  m_SpotJoint = -1;
  numSkeletons = 0;
  for(int skeletonIndex = 0; skeletonIndex < MAX_SKELS; skeletonIndex++)
  {
    m_pSkeleton[skeletonIndex] = NULL;
    m_pMotion[skeletonIndex] = NULL;
  }
}

DisplaySkeleton::~DisplaySkeleton(void)
{
  Reset();
}


//Draws the world coordinate axis
void DisplaySkeleton::DrawSpotJointAxis(void) 
{
//  GLfloat axisLength = 0.5f;
//  glBegin(GL_LINES);
//  // draw x axis in red, y axis in green, z axis in blue 
//  glColor3f(1.0f, 0.2f, 0.2f);
//  glVertex3f(0.0f, 0.0f, 0.0f);
//  glVertex3f(axisLength, 0.0f, 0.0f);
//
//  glColor3f(0.2f, 1.0f, 0.2f);
//  glVertex3f(0.0f, 0.0f, 0.0f);
//  glVertex3f(0.0f, axisLength, 0.0f);
//
//  glColor3f(0.2f, 0.2f, 1.0f);
//  glVertex3f(0.0f, 0.0f, 0.0f);
//  glVertex3f(0.0f, 0.0f, axisLength);
//
//  glEnd();
}

//Build display lists for bones
void DisplaySkeleton::SetDisplayList(int skeletonID, Bone *bone, GLuint *pBoneList)
{
//  GLUquadricObj *qobj;
//  int numbones = m_pSkeleton[skeletonID]->numBonesInSkel(bone[0]);
//  *pBoneList = glGenLists(numbones);
//  qobj=gluNewQuadric();
//
//  gluQuadricDrawStyle(qobj, (GLenum) GLU_FILL);
//  gluQuadricNormals(qobj, (GLenum) GLU_SMOOTH);
//
//  float ambientFskeleton = 0.1f;
//  float diffuseFskeleton = 0.9f;
//  float specularFskeleton = 0.1f;
//
//  int colorIndex = numSkeletons % NUMBER_JOINT_COLORS;
//  float jointShininess = 120.0f;
//  float jointAmbient[4] = {ambientFskeleton * jointColors[colorIndex][0], ambientFskeleton * jointColors[colorIndex][1], ambientFskeleton * jointColors[colorIndex][2], 1.0};
//  float jointDiffuse[4] = {diffuseFskeleton * jointColors[colorIndex][0], diffuseFskeleton * jointColors[colorIndex][1], diffuseFskeleton * jointColors[colorIndex][2], 1.0};
//  float jointSpecular[4] = {specularFskeleton * jointColors[colorIndex][0], specularFskeleton * jointColors[colorIndex][1], specularFskeleton * jointColors[colorIndex][2], 1.0};
//
//  float boneColor[3] = {1.0f, 1.0f, 1.0f};
//  float boneShininess = 120.0f; 
//  float boneAmbient[4] = {ambientFskeleton * boneColor[0], ambientFskeleton * boneColor[1], ambientFskeleton * boneColor[2], 1.0};
//  float boneDiffuse[4] = {diffuseFskeleton * boneColor[0], diffuseFskeleton * boneColor[1], diffuseFskeleton * boneColor[2], 1.0};
//  float boneSpecular[4] = {specularFskeleton * boneColor[0], specularFskeleton * boneColor[1], specularFskeleton * boneColor[2], 1.0};
//
//  double jointRadius = 0.10;
//  double boneRadius = 0.10;
//  double sizeDifferenceJointAndBone = 0.05;
//
//  for(int j=0;j<numbones;j++)
//  {
//    glNewList(*pBoneList + j, GL_COMPILE);
//    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, jointAmbient);
//    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, jointDiffuse);
//    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, jointSpecular);
//    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, jointShininess);
//    glPushMatrix();
//    glScalef(float(bone[j].aspy + sizeDifferenceJointAndBone), float(bone[j].aspy + sizeDifferenceJointAndBone), float(bone[j].aspy + sizeDifferenceJointAndBone));
//    //gluSphere(qobj, jointRadius, 20, 20);
//    glPopMatrix();
//
//    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, boneAmbient);
//    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, boneDiffuse);
//    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, boneSpecular);
//    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, boneShininess);
//    glPushMatrix();
//    glScalef(float(bone[j].aspx), float(bone[j].aspy), 1.0f);
//    //gluCylinder(qobj, boneRadius, boneRadius, bone[j].length, 20, 20);
//
//    // Two disks to close the cylinder at the bottom and the top
//    //gluDisk(qobj, 0.0, boneRadius, 20, 20);
//    glPushMatrix();
//    glTranslatef(0.0f, 0.0f, float(bone[j].length));
//    //gluDisk(qobj, 0.0, boneRadius, 20, 20);
//    glPopMatrix();
//    
//    glPopMatrix();    
//    glEndList();
//  }
}

/*
  Define M_k = Modelview matrix at the kth node (bone) in the heirarchy
  M_k stores the transformation matrix of the kth bone in world coordinates
  Our goal is to draw the (k+1)th bone, using its local information and M_k

  In the k+1th node, compute the following matrices:
  rot_parent_current: this is the rotation matrix that 
  takes us from k+1 to the kth local coordinate system 
  R_k+1 : Rotation matrix for the k+1 th node (bone)
  using angles specified by the AMC file in local coordinates
  T_k+1 : Translation matrix for the k+1th node
  
  The update relation is given by:
  M_k+1 = M_k * (rot_parent_current) * R_k+1 + T_k+1
*/

void DisplaySkeleton::DrawBone(Bone *pBone,int skelNum)
{
  static double z_dir[3] = {0.0, 0.0, 1.0};
  double r_axis[3], theta;

  //Transform (rotate) from the local coordinate system of this bone to it's parent
  //This step corresponds to doing: ModelviewMatrix = M_k * (rot_parent_current)
  glMultMatrixd((double*)&pBone->rot_parent_current);     

  //Draw the local coordinate system for the selected bone.
  if((renderMode == BONES_AND_LOCAL_FRAMES) && (pBone->idx == m_SpotJoint))
  {
    /*
    GLint lightingStatus;
    glGetIntegerv(GL_LIGHTING, &lightingStatus);
    glDisable(GL_LIGHTING);
    DrawSpotJointAxis();
    if (lightingStatus)
      glEnable(GL_LIGHTING);
      */
  }

  //translate AMC (rarely used)
  if(pBone->doftz) 
    glTranslatef(0.0f, 0.0f, float(pBone->tz));
  if(pBone->dofty) 
    glTranslatef(0.0f, float(pBone->ty), 0.0f);
  if(pBone->doftx) 
    glTranslatef(float(pBone->tx), 0.0f, 0.0f);

  //rotate AMC 
  if(pBone->dofrz) 
    glRotatef(float(pBone->rz), 0.0f, 0.0f, 1.0f);
  if(pBone->dofry) 
    glRotatef(float(pBone->ry), 0.0f, 1.0f, 0.0f);
  if(pBone->dofrx) 
    glRotatef(float(pBone->rx), 1.0f, 0.0f, 0.0f);

  //Store the current ModelviewMatrix (before adding the translation part)
  glPushMatrix();

  //Compute tx, ty, tz : translation from pBone to its child (in local coordinate system of pBone)
  double tx = pBone->dir[0] * pBone->length;
  double ty = pBone->dir[1] * pBone->length;
  double tz = pBone->dir[2] * pBone->length;

  // Use the current ModelviewMatrix to display the current bone
  // Rotate the bone from its canonical position (elongated sphere 
  // with its major axis parallel to X axis) to its correct orientation
  if(pBone->idx == Skeleton::getRootIndex())
  {
    // do nothing
  }
  else
  { 
    //Compute the angle between the canonical pose and the correct orientation 
    //(specified in pBone->dir) using cross product.
    //Using the formula: r_axis = z_dir x pBone->dir
    r_axis[0] = z_dir[1]*pBone->dir[2]-z_dir[2]*pBone->dir[1];
    r_axis[1] = z_dir[2]*pBone->dir[0]-z_dir[0]*pBone->dir[2];
    r_axis[2] = z_dir[0]*pBone->dir[1]-z_dir[1]*pBone->dir[0];

    double dot_prod = z_dir[0] * pBone->dir[0] + z_dir[1] * pBone->dir[1] + z_dir[2] * pBone->dir[2] ;
    double r_axis_len = sqrt(r_axis[0] * r_axis[0] + r_axis[1] * r_axis[1] + r_axis[2] * r_axis[2]);
    theta = atan2(r_axis_len, dot_prod);

    glRotatef(float(theta*180./M_PI), float(r_axis[0]), float(r_axis[1]), float(r_axis[2]));
    //glScalef(pBone->aspx, pBone->aspy, 1.0);

    MATRIX4 scaling = MATRIX4::I();
    scaling(0,0) = pBone->aspx;
    scaling(1,1) = pBone->aspy;
    boneScalings[pBone->idx] = scaling;

    float modelview[16];
    glGetFloatv(GL_MODELVIEW, &modelview[0]);
    
    MATRIX4 currentTransform;
    int i = 0;
    for (int x = 0; x < 4; x++)
      for (int y = 0; y < 4; y++, i++)
        currentTransform(x,y) = modelview[i];

    VEC4F currentTranslation;
    currentTranslation[0] = currentTransform(3,0);
    currentTranslation[1] = currentTransform(3,1);
    currentTranslation[2] = currentTransform(3,2);
    currentTransform(3,0) = 0;
    currentTransform(3,1) = 0;
    currentTransform(3,2) = 0;
    //cout << " translation: " << currentTranslation << endl;

    boneRotations[pBone->idx] = currentTransform;
    boneTranslations[pBone->idx] = currentTranslation;
  }

  glPopMatrix(); 

  // Finally, translate the bone, depending on its length and direction
  // This step corresponds to doing: M_k+1 = ModelviewMatrix += T_k+1
  glTranslatef(float(tx), float(ty), float(tz));
}
void DisplaySkeleton::DrawBoneManual(Bone *pBone,int skelNum, MATRIX4& manualTransform)
{
  static double z_dir[3] = {0.0, 0.0, 1.0};
  double r_axis[3], theta;

  //Transform (rotate) from the local coordinate system of this bone to it's parent

  MATRIX4 pBoneRot((double*)&pBone->rot_parent_current);    
  manualTransform = manualTransform * pBoneRot.transpose();


  //Draw the local coordinate system for the selected bone.
  if((renderMode == BONES_AND_LOCAL_FRAMES) && (pBone->idx == m_SpotJoint))
  {
    /*
    GLint lightingStatus;
    glGetIntegerv(GL_LIGHTING, &lightingStatus);
    glDisable(GL_LIGHTING);
    DrawSpotJointAxis();
    if (lightingStatus)
      glEnable(GL_LIGHTING);
      */
  }

  //translate AMC (rarely used)
  if(pBone->doftz){
    manualTransform = manualTransform * MATRIX4::Translation(VEC3F(0.0f, 0.0f, float(pBone->tz)));
  }
  if(pBone->dofty){
    manualTransform = manualTransform * MATRIX4::Translation(VEC3F(0.0f, float(pBone->ty), 0.0f));
  }
  if(pBone->doftx){
    manualTransform = manualTransform * MATRIX4::Translation(VEC3F(float(pBone->tx), 0.0f, 0.0f));
  }

  //rotate AMC 
  if(pBone->dofrz){
    manualTransform = manualTransform * MATRIX4::RotationZ(pBone->rz * M_PI / 180);
  }
  if(pBone->dofry){
    manualTransform = manualTransform * MATRIX4::RotationY(pBone->ry * M_PI / 180);
  }if(pBone->dofrx){
    manualTransform = manualTransform * MATRIX4::RotationX(pBone->rx* M_PI / 180);
  }

  //Store the current ModelviewMatrix (before adding the translation part)
  MATRIX4 matrixCopy(manualTransform);

  //Compute tx, ty, tz : translation from pBone to its child (in local coordinate system of pBone)
  double tx = pBone->dir[0] * pBone->length;
  double ty = pBone->dir[1] * pBone->length;
  double tz = pBone->dir[2] * pBone->length;

  // Use the current ModelviewMatrix to display the current bone
  // Rotate the bone from its canonical position (elongated sphere 
  // with its major axis parallel to X axis) to its correct orientation
  if(pBone->idx == Skeleton::getRootIndex())
  {
    // do nothing
  }
  else
  { 
    //Compute the angle between the canonical pose and the correct orientation 
    //(specified in pBone->dir) using cross product.
    //Using the formula: r_axis = z_dir x pBone->dir
    r_axis[0] = z_dir[1]*pBone->dir[2]-z_dir[2]*pBone->dir[1];
    r_axis[1] = z_dir[2]*pBone->dir[0]-z_dir[0]*pBone->dir[2];
    r_axis[2] = z_dir[0]*pBone->dir[1]-z_dir[1]*pBone->dir[0];

    double dot_prod = z_dir[0] * pBone->dir[0] + z_dir[1] * pBone->dir[1] + z_dir[2] * pBone->dir[2] ;
    double r_axis_len = sqrt(r_axis[0] * r_axis[0] + r_axis[1] * r_axis[1] + r_axis[2] * r_axis[2]);
    theta = atan2(r_axis_len, dot_prod);


    manualTransform = manualTransform * MATRIX4::Rotation(VEC3F(float(r_axis[0]), float(r_axis[1]), float(r_axis[2])), theta);

    MATRIX4 scaling = MATRIX4::I();
    scaling(0,0) = pBone->aspx;
    scaling(1,1) = pBone->aspy;
    boneScalings[pBone->idx] = scaling;

    MATRIX4& currentTransform = manualTransform;

    VEC4F currentTranslation;
    currentTranslation[0] = currentTransform(3,0);
    currentTranslation[1] = currentTransform(3,1);
    currentTranslation[2] = currentTransform(3,2);
    currentTransform(3,0) = 0;
    currentTransform(3,1) = 0;
    currentTransform(3,2) = 0;
    //cout << " translation: " << currentTranslation << endl;

    boneRotations[pBone->idx] = currentTransform;
    boneTranslations[pBone->idx] = currentTranslation;
  }

  manualTransform = matrixCopy * MATRIX4::Translation(VEC3F(float(tx), float(ty), float(tz)));
}
void DisplaySkeleton::SetShadowingModelviewMatrix(double ground[4], double light[4])
{
  double dot;
  double shadowMat[4][4];

  dot = ground[0] * light[0] + ground[1] * light[1] + ground[2] * light[2] + ground[3] * light[3];

  shadowMat[0][0] = dot - light[0] * ground[0];
  shadowMat[1][0] = 0.0 - light[0] * ground[1];
  shadowMat[2][0] = 0.0 - light[0] * ground[2];
  shadowMat[3][0] = 0.0 - light[0] * ground[3];

  shadowMat[0][1] = 0.0 - light[1] * ground[0];
  shadowMat[1][1] = dot - light[1] * ground[1];
  shadowMat[2][1] = 0.0 - light[1] * ground[2];
  shadowMat[3][1] = 0.0 - light[1] * ground[3];

  shadowMat[0][2] = 0.0 - light[2] * ground[0];
  shadowMat[1][2] = 0.0 - light[2] * ground[1];
  shadowMat[2][2] = dot - light[2] * ground[2];
  shadowMat[3][2] = 0.0 - light[2] * ground[3];

  shadowMat[0][3] = 0.0 - light[3] * ground[0];
  shadowMat[1][3] = 0.0 - light[3] * ground[1];
  shadowMat[2][3] = 0.0 - light[3] * ground[2];
  shadowMat[3][3] = dot - light[3] * ground[3];

  glMultMatrixd((const GLdouble*)shadowMat);
}

//Traverse the hierarchy starting from the root 
//Every node in the data structure has just one child pointer. 
//If there are more than one children for any node, they are stored as sibling pointers
//The algorithm draws the current node (bone), visits its child and then visits siblings
void DisplaySkeleton::Traverse(Bone *ptr,int skelNum)
{
  if(ptr != NULL)
  {
    glPushMatrix();
    DrawBone(ptr,skelNum);
    Traverse(ptr->child,skelNum);
    glPopMatrix();
    Traverse(ptr->sibling,skelNum);
  }
}
void DisplaySkeleton::TraverseManual(Bone *ptr, int skelNum, MATRIX4& manualTransform)
{
  if(ptr != NULL)
  {
    MATRIX4 matrixCopy(manualTransform);
    DrawBoneManual(ptr,skelNum, manualTransform);
    TraverseManual(ptr->child,skelNum, manualTransform);
    TraverseManual(ptr->sibling,skelNum, matrixCopy);
  }
}
//Draw the skeleton
void DisplaySkeleton::ComputeBonePositions(RenderMode renderMode_)
{
  unsigned int numbones = m_pSkeleton[0]->numBonesInSkel(*m_pSkeleton[0]->getRoot());
  if (boneRotations.size() != numbones)
  {
    boneRotations.resize(numbones); 
    boneTranslations.resize(numbones); 
    boneScalings.resize(numbones); 
    boneLengths.resize(numbones); 
    for (unsigned int x = 0; x < numbones; x++)
      boneLengths[x] = m_pSkeleton[0]->getBone(x).length;
  }

  // Set render mode
  renderMode = renderMode_;
 
  glPushMatrix();

  // load up identity so we just get the pure transforms out at the end,
  // not the modelview as well
  glLoadIdentity();

  MATRIX4 identity = MATRIX4::I();

  //Translate the root to the correct position (it is (0,0,0) if no motion is loaded)
  //   glTranslatef(m_pSkeleton->m_RootPos[0], m_pSkeleton->m_RootPos[1], m_pSkeleton->m_RootPos[2]);

  //draw the skeleton starting from the root
  for (int i = 0; i < numSkeletons; i++)
  {
    glPushMatrix();
    MATRIX4 currentTransform = identity;

    double translation[3];
    m_pSkeleton[i]->GetTranslation(translation);
    double rotationAngle[3];
    m_pSkeleton[i]->GetRotationAngle(rotationAngle);


    glTranslatef(float(MOCAP_SCALE * translation[0]), float(MOCAP_SCALE * translation[1]), float(MOCAP_SCALE * translation[2]));
    glRotatef(float(rotationAngle[0]), 1.0f, 0.0f, 0.0f);
    glRotatef(float(rotationAngle[1]), 0.0f, 1.0f, 0.0f);
    glRotatef(float(rotationAngle[2]), 0.0f, 0.0f, 1.0f);

    currentTransform = currentTransform * MATRIX4::Translation(MOCAP_SCALE * VEC3F(float(translation[0]), float(translation[1]), float(translation[2])));
    currentTransform = currentTransform * MATRIX4::RotationX(rotationAngle[0] * M_PI / 180);
    currentTransform = currentTransform * MATRIX4::RotationY(rotationAngle[1] * M_PI / 180);
    currentTransform = currentTransform * MATRIX4::RotationZ(rotationAngle[2] * M_PI / 180);

    // Traverse(m_pSkeleton[i]->getRoot(),i);
    TraverseManual(m_pSkeleton[i]->getRoot(),i, currentTransform);

    glPopMatrix();
  }
  glPopMatrix();
}
void DisplaySkeleton::ComputeBonePositionsManual(RenderMode renderMode_)
{
  unsigned int numbones = m_pSkeleton[0]->numBonesInSkel(*m_pSkeleton[0]->getRoot());
  if (boneRotations.size() != numbones)
  {
    boneRotations.resize(numbones); 
    boneTranslations.resize(numbones); 
    boneScalings.resize(numbones); 
    boneLengths.resize(numbones); 
    for (unsigned int x = 0; x < numbones; x++)
      boneLengths[x] = m_pSkeleton[0]->getBone(x).length;
  }

  // Set render mode
  renderMode = renderMode_;

  // load up identity so we just get the pure transforms out at the end,
  // not the modelview as well
  MATRIX4 identity = MATRIX4::I();

  //Translate the root to the correct position (it is (0,0,0) if no motion is loaded)
  //   glTranslatef(m_pSkeleton->m_RootPos[0], m_pSkeleton->m_RootPos[1], m_pSkeleton->m_RootPos[2]);

  //draw the skeleton starting from the root
  for (int i = 0; i < numSkeletons; i++)
  {
    MATRIX4 currentTransform = identity;

    double translation[3];
    m_pSkeleton[i]->GetTranslation(translation);
    double rotationAngle[3];
    m_pSkeleton[i]->GetRotationAngle(rotationAngle);

    currentTransform = currentTransform * MATRIX4::Translation(MOCAP_SCALE * VEC3F(float(translation[0]), float(translation[1]), float(translation[2])));
    currentTransform = currentTransform * MATRIX4::RotationX(rotationAngle[0] * M_PI / 180);
    currentTransform = currentTransform * MATRIX4::RotationY(rotationAngle[1] * M_PI / 180);
    currentTransform = currentTransform * MATRIX4::RotationZ(rotationAngle[2] * M_PI / 180);

    // Traverse(m_pSkeleton[i]->getRoot(),i);
    TraverseManual(m_pSkeleton[i]->getRoot(),i, currentTransform);
  }
}
void DisplaySkeleton::LoadMotion(Motion * pMotion)
{
  // always load the motion for the latest skeleton
  if(m_pMotion[numSkeletons - 1] != NULL) 
    delete m_pMotion[numSkeletons - 1];
  m_pMotion[numSkeletons - 1] = pMotion;
}

//Set skeleton for display
void DisplaySkeleton::LoadSkeleton(Skeleton *pSkeleton)
{
  if (numSkeletons >= MAX_SKELS) 
    return;

  m_pSkeleton[numSkeletons] = pSkeleton;

  //Create the display list for the skeleton
  //All the bones are the elongated spheres centered at (0,0,0).
  //The axis of elongation is the X axis.
  SetDisplayList(numSkeletons, m_pSkeleton[numSkeletons]->getRoot(), &m_BoneList[numSkeletons]);
  numSkeletons++;
}

Motion * DisplaySkeleton::GetSkeletonMotion(int skeletonIndex)
{
  if (skeletonIndex < 0 || skeletonIndex >= MAX_SKELS)
  {
    printf("Error in DisplaySkeleton::GetSkeletonMotion: index %d is illegal.\n", skeletonIndex);
    exit(0);
  }
  return m_pMotion[skeletonIndex];
}

Skeleton * DisplaySkeleton::GetSkeleton(int skeletonIndex)
{
  if (skeletonIndex < 0 || skeletonIndex >= numSkeletons)
  {
    printf("Error in DisplaySkeleton::GetSkeleton: skeleton index %d is illegal.\n", skeletonIndex);
    exit(0);
  }
  return m_pSkeleton[skeletonIndex];
}

void DisplaySkeleton::Reset(void)
{
  for(int skeletonIndex = 0; skeletonIndex < MAX_SKELS; skeletonIndex++)
  {
    if (m_pSkeleton[skeletonIndex] != NULL)
    {
      delete (m_pSkeleton[skeletonIndex]);
      glDeleteLists(m_BoneList[skeletonIndex], 1);
      m_pSkeleton[skeletonIndex] = NULL;
    }
    if (m_pMotion[skeletonIndex] != NULL)
    {
      delete (m_pMotion[skeletonIndex]);
      m_pMotion[skeletonIndex] = NULL;
    }
  }
  numSkeletons = 0;
}



