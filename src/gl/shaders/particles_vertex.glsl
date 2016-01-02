#version 130

uniform mat4 model_view ;
uniform mat4 projection ;

in vec3 vertex ;
in float density ;
in mat4 frame ;

out vec3 normal ;
out vec3 pos ;
out float alpha ;

void main()
{
	vec4 ip = model_view * frame * vec4( vertex, 1)  ;
	pos     = ip.xyz ;
	normal  = normalize( mat3(model_view) * mat3(frame) * vertex ) ;

	alpha = density ;

	gl_Position = projection * ip ;
}
