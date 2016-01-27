#version 130

in vec2 coord ;
in mat2 frame ;
in float val ;

out vec4 color ;

void main (void)
{
    if(val < 1.e-6) discard ;

    float mag = dot( coord, frame*coord );
    if (mag > val-1.e-6) discard;   // kill pixels outside circle

    color = vec4(val, .5, 1-val, 1) ;

}
