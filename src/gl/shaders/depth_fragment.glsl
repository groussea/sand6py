#version 330 core

in vec4 shadow_coord  ;
out float fragmentdepth;

void main(){
//    fragmentdepth = 0.1 ;//shadow_coord.z;
//    gl_FragDepth =  gl_FragCoord.z;
    gl_FragDepth =  shadow_coord.z / 100 ;
//    fragmentdepth =
    //gl_FragColor = vec4(0,fragmentdepth,0,1) ;
}
