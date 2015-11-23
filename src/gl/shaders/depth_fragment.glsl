#version 330 core

in vec4 shadow_coord  ;

void main(){
//    gl_FragDepth = shadow_coord.xy / shadow_coord.w ;
//    gl_FragDepth =  gl_FragCoord.z/10;
    gl_FragDepth =  shadow_coord.z / shadow_coord.w ;
//    fragmentdepth =
    //gl_FragColor = vec4(0,fragmentdepth,0,1) ;
}
