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

in vec3 normal;
in vec3 pos ;
in float alpha ;

out vec4 color ;

void main (void)
{
    vec4 ambientMat = vec4(0.3,0.2, 0.1, 1. );
    vec4 diffuseMat = vec4(0.6, 0.43, 0.38, 1. );
    vec4 specMat    = vec4(1. , 1., 1., 1. );
    float specPow = 15.0;

    vec4 diffuse;
    vec4 spec;
    vec4 ambient;

    vec4 light = vec4( light_pos, 1 )  ;
    vec3 L = normalize( light.xyz - pos);
    vec3 E = normalize(-pos);
    vec3 R = normalize(reflect(-L,normal));

    ambient = ambientMat;
    diffuse = clamp( diffuseMat * max(dot(normal,L), 0.0)  , 0.0, 1.0 ) ;
    spec = clamp ( specMat * pow(max(dot(R,E),0.0),0.3*specPow) , 0.0, 1.0 );

    color = ambient + diffuse + spec;
    color.a = alpha ;

}
