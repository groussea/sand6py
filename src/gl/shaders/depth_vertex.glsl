#version 130

uniform mat4 depth_mvp;

in vec3 vertex;

out vec4 shadow_coord  ;

void main(){
  gl_Position =  depth_mvp * vec4(vertex,1);
  shadow_coord = gl_Position ;
}
