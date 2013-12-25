Peter Bennion & Owen Campbell
CS180: Computer Graphics
Project 4: Diffuse Raytracing

Submission is designed to be run on the computers at
UCSB's CSIL tech lab. Makefile and includes are for a
Linux Fedora distribution, with GLU and GLUT installed.

This is the final assignment for CS180, built upon my
earlier raytracer.

MATRIX3.h, MATRIX4.h, VEC3F.h, VEC4F.h, types.h,
skeleton.h, posture.h, motion.h, displaySkeleton.h,
bonevector.h, and their corresponding .cpp's were
provided by the professor as a framework for
animation.

The purpose of this assignment was to add distributed
raytracing effects -- such as soft shadows, imperfect 
reflections, and textures. The former two effects use
a stratified sampling approach to acquire random samples
for the shadow and reflection rays.

A large portion of the new code was added by my partner,
Owen Campbell.