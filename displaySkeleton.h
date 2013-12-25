/*
display.h

Display the skeleton, ground plane and other objects.			

Revision 1 - Steve Lin, Jan. 14, 2002
Revision 2 - Alla and Kiran, Jan 18, 2002
Revision 3 - Jernej Barbic and Yili Zhao, Feb, 2012
*/

#ifndef _DISPLAY_SKELETON_H_
#define _DISPLAY_SKELETON_H_

//osx
//#include <OpenGL/glu.h>

//linux
#include <GL/glu.h>


#include "skeleton.h"
#include "motion.h"
#include <vector>
#include <map>
#include "MATRIX4.h"

using namespace std;

class DisplaySkeleton 
{

//member functions
public: 
  enum RenderMode
  {
    BONES_ONLY, BONES_AND_LOCAL_FRAMES
  };
  enum JointColor
  {
    GREEN, RED, BLUE, NUMBER_JOINT_COLORS
  };

  DisplaySkeleton();
  ~DisplaySkeleton();

  //set skeleton for display
  void LoadSkeleton(Skeleton * pSkeleton);
  //set motion for display
  void LoadMotion(Motion * pMotion);

  //display the scene (skeleton, ground plane ....)
  void ComputeBonePositions(RenderMode renderMode);
  void ComputeBonePositionsManual(RenderMode renderMode);

  void SetDisplayedSpotJoint(int jointID) {m_SpotJoint = jointID;}
  int GetDisplayedSpotJoint(void) {return m_SpotJoint;}
  int GetNumSkeletons(void) {return numSkeletons;}
  Skeleton * GetSkeleton(int skeletonIndex);
  Motion * GetSkeletonMotion(int skeletonIndex);

  void Reset(void);

  vector<MATRIX4>& rotations()  { return boneRotations; };
  vector<MATRIX4>& scalings()   { return boneScalings; };
  vector<VEC4F>& translations() { return boneTranslations; };
  vector<float>& lengths()      { return boneLengths; };

protected:
  RenderMode renderMode;
  // Draw a particular bone
  void DrawBone(Bone *ptr, int skelNum);
  // Draw the skeleton hierarchy
  void Traverse(Bone *ptr, int skelNum);

  // Draw a particular bone
  void DrawBoneManual(Bone *ptr, int skelNum, MATRIX4& manualTransform);
  // Draw the skeleton hierarchy
  void TraverseManual(Bone *ptr, int skelNum, MATRIX4& manualTransform);

  // Model matrix for the shadow
  void SetShadowingModelviewMatrix(double ground[4], double light[4]);
  void DrawSpotJointAxis(void);
  void SetDisplayList(int skeletonID, Bone *bone, GLuint *pBoneList);

  int m_SpotJoint;		//joint whose local coordinate framework is drawn
  int numSkeletons;
  Skeleton *m_pSkeleton[MAX_SKELS];		//pointer to current skeleton
  Motion *m_pMotion[MAX_SKELS];		//pointer to current motion	
  GLuint m_BoneList[MAX_SKELS];		//display list with bones

  static float jointColors[NUMBER_JOINT_COLORS][3];

  //std::map<Bone*, MATRIX4> boneTransforms;
  vector<MATRIX4> boneRotations;
  vector<MATRIX4> boneScalings;
  vector<VEC4F> boneTranslations;
  vector<float> boneLengths;
};

#endif
