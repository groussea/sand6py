#version 130

in vec3 vertex;

uniform mat4 depth_mvp;

out vec4 shadow_coord  ;

void main(){

 gl_Position =  depth_mvp * vec4(vertex,1);
 shadow_coord = gl_Position ;
}
