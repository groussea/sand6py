#version 130

uniform mat4 model_view ;
uniform mat4 projection ;
uniform float radius ;
uniform vec3 center ;

in vec3 vertex ;

out vec2 coord ;

void main()
{
    coord = vertex.xy ;
    vec4 depl = transpose(model_view) * vec4(vertex, 1) ;

    vec3 pos = center + depl.xyz*radius ;


    gl_Position = projection * model_view * vec4(pos,1) ;
}


