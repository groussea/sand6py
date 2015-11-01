#version 140

uniform mat4 model_view ;
uniform mat4 projection ;

in vec3 vertex ;
//varying vec3 pos ;

out vec3 normal ;
out vec3 pos ;

void main()
{
//  pos = vec3(gl_ModelViewMatrix * gl_Vertex);
//  normal = normalize( gl_NormalMatrix * gl_Vertex.xyz );


	normal = vertex ;
	vec3 ip = vertex + gl_InstanceID * vec3( 1., 0, 0 ) ;
	pos = ( model_view * vec4( ip, 1 ) ).xyz ;

	gl_Position = projection * vec4( pos, 1) ;
}
