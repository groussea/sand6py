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

uniform mat4 depth_mvp;
uniform float grain_size  ;
uniform float pixel_size  ;

in vec3 vertex;
in float visibility ;

out vec4 shadow_coord  ;
out float vis ;

void main(){
 vis = visibility ;

 gl_Position =  depth_mvp * vec4(vertex,1);
 shadow_coord = gl_Position ;

 gl_PointSize = max(1, pixel_size * grain_size / gl_Position.w );
}
