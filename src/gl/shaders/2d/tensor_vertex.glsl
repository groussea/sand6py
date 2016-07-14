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
in vec3 value ;

out vec2 coord ;
out mat2 frame ;
out float val ;

void main()
{
    val = length(value) ;

    coord.y = 1 - 2*((gl_VertexID>>1)%2) ;
    coord.x = 1 - 2*(((gl_VertexID>>1)%2) ^ (gl_VertexID%2)) ;

    vec2 vertex_world = vertex + val * coord ;

    gl_Position = projection * model_view * vec4( vertex_world, 0, 1) ;

    coord = coord ;

    frame = 1/( (value.x + value.y) * (value.x - value.y) - value.z*value.z )
            * mat2( vec2(value.x-value.y, -value.z), vec2(-value.z, value.x+value.y) ) ;
}
