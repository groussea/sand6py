#version 130

uniform mat4 model_view ;
uniform mat4 projection ;

in vec2 vertex ;
in vec3 value ;

out vec2 coord ;
out mat2 frame ;
out float val ;

void main()
{
    val = length(value) ;

    coord.y = 1 - 2*((gl_VertexID>>1)%2) ;
    coord.x = 1 - 2*(((gl_VertexID>>1)%2) ^ (gl_VertexID%2)) ;

    vec2 vertex_world = vertex + val * coord ;

    gl_Position = projection * model_view * vec4( vertex_world, 0, 1) ;

    coord = coord ;

    frame = 1/( (value.x + value.y) * (value.x - value.y) - value.z*value.z )
            * mat2( vec2(value.x-value.y, -value.z), vec2(-value.z, value.x+value.y) ) ;
}
