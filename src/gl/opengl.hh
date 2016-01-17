
#ifndef GL_GLEXT_PROTOTYPES
#define GL_GLEXT_PROTOTYPES
#endif

#ifdef __APPLE__
#include <OpenGL/gl.h>
#define  glVertexAttribDivisor glVertexAttribDivisorARB
#define  glDrawArraysInstanced glDrawArraysInstancedARB
#define  GL_PROGRAM_POINT_SIZE GL_PROGRAM_POINT_SIZE_EXT
#else
#include <GL/gl.h>
#endif
