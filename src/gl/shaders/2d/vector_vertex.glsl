#version 130

uniform mat4 model_view ;
uniform mat4 projection ;

in vec2 vertex ;
in vec2 value ;

out float val ;

void main()
{
    vec3 bary = vec3(0,0,0) ;
    bary[gl_VertexID] = 1. ;

    vec2 n = .25*vec2( -value.y, value.x ) ;

    vec2 P1 = .66*value ;
    vec2 P2 = -n -.33*value ;
    vec2 P3 =  n -.33*value ;

    vec2 vertex_world = vertex + bary.x*P1 + bary.y*P2 + bary.z*P3 ;

    gl_Position = projection * model_view * vec4( vertex_world, 0, 1) ;

    val = length(value) ;
}
