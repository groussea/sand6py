#version 130

uniform mat4 model_view ;
uniform mat3 rotation ;
uniform float radius ;
uniform vec3 center ;
uniform vec3 light_pos ;

in vec2 coord ;

out vec4 color ;

void main (void)
{

    float mag = dot( coord, coord );
    if (mag > 1.0) discard;   // kill pixels outside circle

    vec3 N = vec3( coord, sqrt(1 - mag) ) ;

    vec3 normal_w = transpose(rotation) * mat3( transpose(model_view) ) * N ;
    vec3 pos_w = center + transpose(rotation) *N*radius ;

    vec3 pos = ( model_view * vec4( pos_w, 1)).xyz ;
    vec3 normal = normalize( mat3(model_view) * rotation * normal_w );

//    color = vec4(coord,0,1)   ;
	
	vec4 ambientMat = vec4(0.4 ,  0.1, 0.1, 1. );
	vec4 diffuseMat = vec4(0.45, 0.15, 0.1, 1. );
	vec4 specMat    = vec4(0.5 ,  0.3, 0.1, 1. );
	
	if( normal_w.x > 0 ) {
		ambientMat = vec4(0.3 , 0.1,  0.1, 1. );
		diffuseMat = vec4(0.35, 0.1, 0.15, 1. );
		specMat    = vec4(0.4 , 0.1,  0.3, 1. );
	}
	if( normal_w.y > 0 ) {
		ambientMat = vec4( 0.1, 0.1, 0.5 , 1. );
		diffuseMat = vec4(0.15, 0.1, 0.55, 1. );
		specMat    = vec4( 0.3, 0.1, 0.7 , 1. );
	}
	if( normal_w.z > 0 ) {
		ambientMat = vec4( 0.1,  0.3, 0.1, 1. );
		diffuseMat = vec4(0.15, 0.35, 0.1, 1. );
		specMat    = vec4(0.3 ,  0.4, 0.1, 1. );
	}


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

}
