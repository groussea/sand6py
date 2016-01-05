#version 130

uniform mat4 depth_mvp;
uniform float grain_size  ;
uniform float pixel_size  ;

in vec3 vertex;
in float visibility ;

out vec4 shadow_coord  ;
out float vis ;

void main(){
 vis = visibility ;

 gl_Position =  depth_mvp * vec4(vertex,1);
 shadow_coord = gl_Position ;

 gl_PointSize = max(1, pixel_size * grain_size / gl_Position.w );
}
