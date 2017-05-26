/*
 * This file is part of Sand6, a C++ continuum-based granular simulator.
 *
 * Copyright 2016 Gilles Daviet <gilles.daviet@inria.fr> (Inria - Université Grenoble Alpes)
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

#version 130

uniform mat4 model_view ;
uniform mat4 projection ;

in vec2 vertex ;
in vec2 value ;

out float val ;

void main()
{
    vec3 bary = vec3(0,0,0) ;
    bary[gl_VertexID] = 1. ;

    vec2 n = .25*vec2( -value.y, value.x ) ;

    vec2 P1 = .66*value ;
    vec2 P2 = -n -.33*value ;
    vec2 P3 =  n -.33*value ;

    vec2 vertex_world = vertex + bary.x*P1 + bary.y*P2 + bary.z*P3 ;

    gl_Position = projection * model_view * vec4( vertex_world, 0, 1) ;

    val = length(value) ;
}
