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

in vec2 coord ;
in float val ;
in float alpha ;

out vec4 color ;

void main (void)
{

	float mag = dot( coord, coord );
	if (mag > 1.0) discard;   // kill pixels outside circle

	vec4 ambientMat  = vec4(0,0,0,0) ;

	float mat_0 = clamp( 1 - 2*abs(val/5-1.0), 0.0, 1.0 ) ;
	float mat_1 = clamp( 1 - 2*abs(val/5-0.5), 0.0, 1.0 ) ;
	float mat_2 = clamp( 1 - 2*abs(val/5-0.0), 0.0, 1.0 ) ;

	float mat_3 = clamp( 1 + 2*   (val/5-1.5), 0.0, 1.0 ) ;

	ambientMat += mat_0*vec4( vec3(1.0, 0.0, 0.0), 1. );
	ambientMat += mat_1*vec4( vec3(0.0, 1.0, 0.0), 1. );
	ambientMat += mat_2*vec4( vec3(0.0, 0.0, 1.0), 1. );
	ambientMat += mat_3*vec4( vec3(1.0, 1.0, 1.0), 1. );

	color = ambientMat ;
	color.a = alpha ;

}
