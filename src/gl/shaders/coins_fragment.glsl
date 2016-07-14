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
uniform vec3 light_pos ;
uniform sampler2D depth_texture;

in vec2 coord ;
in vec3 normal_eye ;
in vec3 vertex_eye ;
in vec4 shadow_coord  ;

in float vis ;
in float material ;

out vec4 color ;

void main (void)
{
	if( vis < 0 || dot(coord.xy, coord.xy) > 1 )
		discard ;

	float alpha = 1 ; //vis ;

	vec2 tex_coords = 0.5 * ( shadow_coord.xy / shadow_coord.w ) + vec2(0.5,0.5) ;

	float zz = shadow_coord.z / shadow_coord.w  ;
	float zs = texture( depth_texture, tex_coords ).r ;

	if ( zs  <  zz ){
		 alpha = 1. - (zz-zs)*25 ;
	 }
	alpha = pow( clamp(alpha,0,1), 3) ;

	vec4 ambientMat  = vec4(0,0,0,0) ;
	vec4 diffuseMat  = vec4(0,0,0,0) ;

	if( material < .4 )
		ambientMat = vec4( 212,175,55, 255.)/255;
	else if( material < .7 )
		ambientMat = vec4( 192,192,192, 255.)/255;
	else
		ambientMat = vec4( 200,117,51, 255.)/255;

	diffuseMat = ambientMat ;
	vec4 specMat    = vec4( .8, .8, .8, 1. );


	float specPow = 15.0;

	vec4 diffuse;
	vec4 spec;
	vec4 ambient;

	vec4 light = model_view * vec4( light_pos, 1 )  ;
	vec3 L = normalize( light.xyz - vertex_eye);
	vec3 E = normalize(- vertex_eye);
	vec3 R = normalize(reflect(-L, normal_eye ));

	ambient = ambientMat;
	diffuse = diffuseMat * clamp( abs(dot(normal_eye,L)) , 0.0, 1.0 ) ;
	spec =  specMat * clamp (pow(abs(dot(R,E)),0.3*specPow) , 0.0, 1.0 );

	color = ( 0.4 + 0.6*alpha) * .5 * ambient ;
	color += ( 0.2 + 0.8*alpha ) * .5 * diffuse ;
	color += alpha*alpha*.5*spec;

	color.a = 1. ;
}
