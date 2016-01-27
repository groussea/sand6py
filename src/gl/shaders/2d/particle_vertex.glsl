#version 130

uniform mat4 model_view ;
uniform mat4 projection ;

in vec2 vertex ;
in float density ;
in float value ;
in vec4 frame ;

out vec2 coord ;
out float val ;
out float alpha ;

void main()
{
    coord.y = 1 - 2*((gl_VertexID>>1)%2) ;
    coord.x = 1 - 2*(((gl_VertexID>>1)%2) ^ (gl_VertexID%2)) ;

    vec2 fcoord = vec2( frame.x * coord.x + frame.z*coord.y, frame.y * coord.x + frame.w*coord.y ) ;
    //mat2 fr = mat2( frame.xy, frame.yz ) ;
    vec2 vertex_world = vertex + fcoord ;

    gl_Position = projection * model_view * vec4( vertex_world, 0, 1) ;

    alpha = density ;
    val = value ;
}
