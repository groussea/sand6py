
/*
 * This file is part of Sand6, a C++ continuum-based granular simulator.
 *
 * Copyright 2016 Gilles Daviet <gilles.daviet@inria.fr> (Inria - Universit√© Grenoble Alpes)
 *
 * Sand6 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * Sand6 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with Sand6.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GL_GLEXT_PROTOTYPES
#define GL_GLEXT_PROTOTYPES
#endif

#ifdef __APPLE__

#define GL_CORE
#ifdef GL_CORE
#include <OpenGL/gl3.h>
#else
#include <OpenGL/gl.h>
#define  glVertexAttribDivisor glVertexAttribDivisorARB
#define  glDrawArraysInstanced glDrawArraysInstancedARB
#define  glDrawElementsInstanced glDrawElementsInstancedARB
#define  GL_PROGRAM_POINT_SIZE GL_PROGRAM_POINT_SIZE_EXT
#endif
#else
#include <GL/gl.h>
#endif
