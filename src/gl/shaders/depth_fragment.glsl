#version 130 

in vec4 shadow_coord  ;

void main(){
//    gl_FragDepth =  gl_FragCoord.z;
    gl_FragDepth =  shadow_coord.z / shadow_coord.w ;
}
