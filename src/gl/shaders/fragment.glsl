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
uniform vec3 ambient ;
uniform sampler2D depth_texture;

in vec3 normal_eye;
in vec3 vertex_eye ;
in vec3 tex_coord ;
in vec4 shadow_coord  ;

out vec4 color ;

void main (void)
{

    //shadow
    float alpha = 1 ; //vis ;

    vec2 tex_coords = 0.5 * ( shadow_coord.xy / shadow_coord.w ) + vec2(0.5,0.5) ;

    float zz = shadow_coord.z / shadow_coord.w  ;
    float zs = texture( depth_texture, tex_coords ).r ;

    if ( zs  <  zz ){
         alpha = 1. - (zz-zs)*10 ;
     }

    alpha = pow( clamp(alpha,0,1), 3) ;

	//other

    vec4 ambientMat = vec4(ambient, 1. );
    vec4 diffuseMat = vec4(.9 , .6, .4, 1. );
    vec4 specMat    = vec4(1. , .8, .7, 1. );

    if( bool( (int(tex_coord.x*10)&1) ^ (int(tex_coord.y*10)&1) ) ) {
        ambientMat *= .7 ;
        diffuseMat *= .7 ;
    }

    float specPow = 12.0;

    vec4 diffuse;
    vec4 spec;

    vec4 light = model_view * vec4( light_pos, 1 )  ;
    vec3 L = normalize( light.xyz - vertex_eye);
    vec3 E = normalize(-vertex_eye);
    vec3 R = normalize(reflect(-L,normal_eye));

    diffuse = clamp( (diffuseMat-ambientMat) * max(dot(normal_eye,L), 0.0)  , 0.0, 1.0 ) ;
    spec = clamp ( (specMat-diffuseMat) * pow(max(dot(R,E),0.0),0.3*specPow) , 0.0, 1.0 );

    color = ambientMat + alpha* diffuse + alpha * alpha * spec;
    color.a = 1 ;
}
