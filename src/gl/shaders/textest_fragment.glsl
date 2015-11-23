#version 130

uniform sampler2D in_texture;

in  vec2 coord ;
out vec4 color ;

void main (void)
{
    vec2 tex_coord = 0.5 * coord + vec2(0.5,0.5) ;
    color = vec4( texture( in_texture, tex_coord ).rgb, 1 ) ;
    //color = vec4( tex_coord, 0, 1 );
}
