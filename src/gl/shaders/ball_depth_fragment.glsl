#version 130

in vec2 coord ;
in vec4 shadow_coord ;

out vec4 color ;

void main (void)
{
    float mag = dot( coord, coord );
    if (mag > 1.0) discard;   // kill pixels outside circle
	
	gl_FragDepth =  shadow_coord.z / shadow_coord.w ;
}
