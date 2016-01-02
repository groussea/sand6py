#version 130

uniform mat4 model_view ;
uniform mat4 projection ;

in vec3 vertex ;
in vec3 normal ;

out vec3 normal_eye ;
out vec3 vertex_eye ;

void main()
{
	vec4 ip = model_view * vec4( vertex, 1)  ;
	vertex_eye  = ip.xyz ;
	normal_eye = normalize( model_view*vec4(normal,0) ).xyz ;

	gl_Position = projection * ip ;
}
