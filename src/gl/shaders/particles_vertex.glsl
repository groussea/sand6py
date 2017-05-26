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

in vec3 vertex ;
in float density ;
in mat4 frame ;

out vec3 normal ;
out vec3 pos ;
out float alpha ;

void main()
{
	vec4 ip = model_view * frame * vec4( vertex, 1)  ;
	pos     = ip.xyz ;
	normal  = normalize( mat3(model_view) * mat3(frame) * vertex ) ;

	alpha = density ;

	gl_Position = projection * ip ;
}
