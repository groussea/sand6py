#version 130

uniform mat4 model_view ;
uniform mat4 projection ;
uniform float radius ;
uniform vec3 center ;

in vec3 vertex ;

out vec2 coord ;
out vec4 shadow_coord ;

void main()
{
    coord = vertex.xy ;
    vec3 depl = transpose(mat3( model_view)) * vertex ;

    vec3 pos = center + depl * radius ;

    gl_Position = projection * model_view * vec4(pos,1) ;
	shadow_coord = gl_Position ;
}


