#version 130

in vec4 shadow_coord  ;

in float vis ;

void main(){
    if( vis < 0 ) discard ;
//    gl_FragDepth =  gl_FragCoord.z;
    gl_FragDepth =  shadow_coord.z / shadow_coord.w ;
}
