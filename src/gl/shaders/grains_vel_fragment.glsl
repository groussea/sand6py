
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


float colormap_green(float x) {
	if (x < 0.09825118520770205) {
		return 5.07556923076926E+02 * x + 1.64923076923077E+02;
	} else if (x < 0.2009111350471108) {
		return 2.86362637362647E+02 * x + 1.86655677655676E+02;
	} else if (x < 0.2994418666360456) {
		return 8.90415982485030E+01 * x + 2.26299671592774E+02;
	} else if (x < 0.5001300871372223) {
		return 9.81627851140242E+00 * x + 2.50023049219689E+02;
	} else if (x < 0.9039205014705658) {
		return ((-3.30848798119696E+01 * x - 5.65722561191396E+02) * x + 2.78046782759626E+02) * x + 2.61515979057614E+02;
	} else {
		return -2.53583846153761E+02 * x + 2.55396153846073E+02;
	}
}

float colormap_red(float x) {
	if (x < 0.1105575469849737) {
		return 4.79433455433456E+02 * x + 3.65079365079361E-01;
	} else if (x < 0.3151890079472769) {
		return 6.25896582484846E+02 * x - 1.58275246854709E+01;
	} else if (x < 0.4023888287265409) {
		return 4.80700000000005E+02 * x + 2.99368421052611E+01;
	} else if (x < 0.5007980763912201) {
		return 3.22042124542111E+02 * x + 9.37789987790044E+01;
	} else if (x < 0.9266376793384552) {
		return ((-2.91150627193739E+02 * x + 2.73891595228739E+02) * x - 1.97954551648389E+02) * x + 3.22069054828072E+02;
	} else {
		return -4.70385384615211E+02 * x + 5.78034615384465E+02;
	}
}

float colormap_blue(float x) {
	if (x < 0.1007720845701718) {
		return 1.66813186813184E+01 * x + 3.72910052910053E+01;
	} else if (x < 0.2891807195246389) {
		return 2.86155895159223E+02 * x + 1.01354904806627E+01;
	} else if (x < 0.4061884871072265) {
		return 4.02182758620675E+02 * x - 2.34172413793071E+01;
	} else if (x < 0.5018816861329155) {
		return 5.35500000000025E+02 * x - 7.75691699604942E+01;
	} else if (x < 0.604070194492165) {
		return -5.10170329670400E+02 * x + 4.47233618233660E+02;
	} else if (x < 0.7060918916718424) {
		return -3.26878215654109E+02 * x + 3.36512315270959E+02;
	} else if (x < 0.812819402403008) {
		return -6.62557264957455E+01 * x + 1.52488888888906E+02;
	} else {
		return -2.16444081632622E+02 * x + 2.74564897959153E+02;
	}
}

vec4 colormap(float x) {
	float r = clamp(colormap_red(x) / 255.0, 0.0, 1.0);
	float g = clamp(colormap_green(x) / 255.0, 0.0, 1.0);
	float b = clamp(colormap_blue(x) / 255.0, 0.0, 1.0);
	return vec4(r, g, b, 1.0);
}


void main (void)
{
	if( vis < 0 ) discard ;

	float alpha = 1 ; //vis ;

	vec4 ambientMat  = vec4(0,0,0,0) ;
	float factor = 2;
	float mat_0 = clamp(1 - 2 * abs(vis / factor  - 1.0), 0.0, 1.0);
	float mat_1 = clamp( 1 - 2*abs(vis/factor -0.5), 0.0, 1.0 ) ;
	float mat_2 = clamp( 1 - 2*abs(vis/factor -0.0), 0.0, 1.0 ) ;

	ambientMat += mat_0*vec4( vec3(0.8, 0.0, 0.0), 1. );
	ambientMat += mat_1*vec4( vec3(0.0, 0.8, 0.0), 1. );
	ambientMat += mat_2*vec4( vec3(0.0, 0.0, 0.8), 1. );

	ambientMat = colormap(vis / 3);
	ambientMat *= 1 + 0.25 * (material - 0.5);
	//    diffuseMat = 1.1 * ambientMat;

	color = ambientMat ;
	// color =
	color.a = 1. ;
}
