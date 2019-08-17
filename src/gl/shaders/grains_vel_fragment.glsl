
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
uniform vec3 light_pos ;
uniform sampler2D depth_texture;

in vec3 normal_eye ;
in vec3 vertex_eye ;
in vec4 shadow_coord  ;

in float vis ;
in float material ;

out vec4 color ;

void main (void)
{
	if( vis < 0 ) discard ;

	float alpha = 1 ; //vis ;

	vec4 ambientMat  = vec4(0,0,0,0) ;

	float mat_0 = clamp( 1 - 2*abs(vis/5-1.0), 0.0, 1.0 ) ;
	float mat_1 = clamp( 1 - 2*abs(vis/5-0.5), 0.0, 1.0 ) ;
	float mat_2 = clamp( 1 - 2*abs(vis/5-0.0), 0.0, 1.0 ) ;

	ambientMat += mat_0*vec4( vec3(1.0, 0.0, 0.0), 1. );
	ambientMat += mat_1*vec4( vec3(0.0, 1.0, 0.0), 1. );
	ambientMat += mat_2*vec4( vec3(0.0, 0.0, 1.0), 1. );

	ambientMat *= 1 + 0.25*(material - 0.5) ;
//    diffuseMat = 1.1 * ambientMat;

	color = ambientMat ;
	color.a = 1. ;
}
