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

#version 330

uniform mat4 model_view ;
uniform mat4 projection ;
uniform mat4 depth_mvp ;

in vec3 vertex ;
in vec3 normal ;
in vec3 uv ;

out vec3 normal_eye ;
out vec3 vertex_eye ;
out vec3 tex_coord ;
out vec4 shadow_coord  ;

void main()
{
    vec4 ip = model_view * vec4( vertex, 1)  ;
    vertex_eye  = ip.xyz ;
    normal_eye = normalize( model_view*vec4(normal,0) ).xyz ;

    gl_Position = projection * ip ;

    tex_coord = uv ;
	shadow_coord = ( depth_mvp * vec4( vertex, 1 ) ) ;
}
