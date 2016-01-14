#version 130

uniform mat4 model_view ;

varying vec3 normal;
varying vec3 pos ;
varying float alpha ;

out vec4 color ;

void main (void)
{
    vec3 light_pos = vec3(0,0,100) ;

    vec4 ambientMat = vec4(0.3,0.2, 0.1, 1. );
    vec4 diffuseMat = vec4(.7 , .6, .2, 1. );
    vec4 specMat    = vec4(1. , 1., 1., 1. );
    float specPow = 15.0;

    vec4 diffuse;
    vec4 spec;
    vec4 ambient;

    vec4 light = model_view * vec4( light_pos, 1 )  ;
    vec3 L = normalize( light.xyz - pos);
    vec3 E = normalize(-pos);
    vec3 R = normalize(reflect(-L,normal));

    ambient = ambientMat;
    diffuse = clamp( diffuseMat * max(dot(normal,L), 0.0)  , 0.0, 1.0 ) ;
    spec = clamp ( specMat * pow(max(dot(R,E),0.0),0.3*specPow) , 0.0, 1.0 );

    color = ambient + diffuse + spec;
    color.a = alpha ;
}
