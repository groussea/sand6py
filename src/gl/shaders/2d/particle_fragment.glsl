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

	float vis = log(1+100*val) ;

	float mat_0 = clamp( 1 - 2*abs(vis/5-1.0), 0.0, 1.0 ) ;
	float mat_1 = clamp( 1 - 2*abs(vis/5-0.5), 0.0, 1.0 ) ;
	float mat_2 = clamp( 1 - 2*abs(vis/5-0.0), 0.0, 1.0 ) ;

	float mat_3 = clamp( 1 + 2*   (vis/5-1.5), 0.0, 1.0 ) ;

	ambientMat += mat_0*vec4( vec3(1.0, 0.0, 0.0), 1. );
	ambientMat += mat_1*vec4( vec3(0.0, 1.0, 0.0), 1. );
	ambientMat += mat_2*vec4( vec3(0.0, 0.0, 1.0), 1. );
	ambientMat += mat_3*vec4( vec3(1.0, 1.0, 1.0), 1. );

	color = ambientMat ;
	color.a = alpha ;

}
