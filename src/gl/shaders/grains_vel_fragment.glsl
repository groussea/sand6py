
#version 130

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
