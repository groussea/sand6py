#version 130

in float val ;

out vec4 color ;

void main (void)
{

//    color = vec4(1,0,0,1) ;
	color = vec4(val,val,val,0.5) ;

}
