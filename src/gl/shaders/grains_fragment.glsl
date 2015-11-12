#version 130

uniform mat4 model_view ;
uniform vec3 light_pos ;

varying vec3 normal_screen ;
varying vec3 vertex_screen ;

varying float alpha ;

out vec4 color ;

void main (void)
{
    vec4 ambientMat = vec4( vec3(0.3,0.3, 0.), 1. );
    vec4 diffuseMat = vec4( alpha*vec3(.8 , .8, 0.), 1. );
    vec4 specMat    = vec4( alpha*alpha*vec3(1. , 1., 1.), 1. );
    float specPow = 15.0;

    vec4 diffuse;
    vec4 spec;
    vec4 ambient;

    vec4 light = model_view * vec4( light_pos, 1 )  ;
    vec3 L = normalize( light.xyz - vertex_screen);
    vec3 E = normalize(- vertex_screen);
    vec3 R = normalize(reflect(-L, normal_screen ));

    ambient = ambientMat;
    diffuse = clamp( diffuseMat * max(dot(normal_screen,L), 0.0)  , 0.0, 1.0 ) ;
    spec = clamp ( specMat * pow(max(dot(R,E),0.0),0.3*specPow) , 0.0, 1.0 );

    color = ambient + diffuse + spec;
//    color.a = alpha ;
}
