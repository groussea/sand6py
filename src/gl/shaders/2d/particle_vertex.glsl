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
in float density ;
in float value ;
in vec4 frame ;

out vec2 coord ;
out float val ;
out float alpha ;

void main()
{
    coord.y = 1 - 2*((gl_VertexID>>1)%2) ;
    coord.x = 1 - 2*(((gl_VertexID>>1)%2) ^ (gl_VertexID%2)) ;

    vec2 fcoord = vec2( frame.x * coord.x + frame.z*coord.y, frame.y * coord.x + frame.w*coord.y ) ;
    //mat2 fr = mat2( frame.xy, frame.yz ) ;
    vec2 vertex_world = vertex + fcoord ;

    gl_Position = projection * model_view * vec4( vertex_world, 0, 1) ;

    alpha = density ;
    val = value ;
}
