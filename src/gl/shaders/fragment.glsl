#version 130

uniform mat4 model_view ;
uniform vec3 light_pos ;
uniform vec3 ambient ;


in vec3 normal_eye;
in vec3 vertex_eye ;
in vec3 tex_coord ;

out vec4 color ;

void main (void)
{
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

    color = ambientMat + diffuse + spec;
    color.a = 1 ;
}
