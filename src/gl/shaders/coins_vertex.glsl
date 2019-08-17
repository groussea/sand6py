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
uniform mat4 depth_mvp  ;
uniform float grain_size  ;

in vec3 vertex ;
in vec3 normal ;
in float visibility ;
in float noise ;

out vec2 coord ;
out vec3 normal_eye ;
out vec3 vertex_eye ;
out vec4 shadow_coord  ;
out float vis ;
out float material ;

void main()
{
    const vec2[6] lut = vec2[6](vec2(-1,-1), vec2(-1,1), vec2(1, 1), vec2(-1, -1), vec2(1, 1), vec2(1, -1) );
    coord = lut[gl_VertexID];

    vec3 t1, t2 ;
    if( normal.x > 1.e-12 ) {
        t1 = normalize( vec3( -normal.y, normal.x, 0 ) ) ;
    } else {
        t1 = normalize( vec3( 0, -normal.z, normal.y ) ) ;
    }
    t2 = cross( normal, t1 );

    vec3 vertex_world = vertex + .5 * grain_size * ( coord.x*t1 + coord.y*t2 ) ;

    vertex_eye = ( model_view * vec4(vertex_world,1) ).xyz ;
    normal_eye  = normalize( model_view * vec4(normal,0) ).xyz ;

    gl_Position = projection * vec4(vertex_eye,1) ;
    shadow_coord = ( depth_mvp * vec4( vertex_world, 1 ) ) ;

    vis = visibility ;
    material = noise ;
}


