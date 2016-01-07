#version 130

in vec2 coord ;
in vec4 shadow_coord  ;
in float vis ;

out vec4 color ;

void main (void)
{
    if( vis < 0 || dot(coord.xy, coord.xy) > 1 )
        discard ;

    gl_FragDepth =  shadow_coord.z / shadow_coord.w ;
}
