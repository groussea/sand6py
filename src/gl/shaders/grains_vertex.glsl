#version 130

uniform mat4 model_view ;
uniform mat4 projection ;

in vec3 vertex ;
in vec3 normal ;

out vec3 normal_screen ;
out vec3 vertex_screen ;
out float alpha ;

void main()
{
	vec4 ip = model_view * vec4( vertex, 1)  ;
	vertex_screen  = ip.xyz ;
	normal_screen  = normalize( mat3(model_view) * normal ) ;

	alpha = 1.0 ;

	gl_Position = projection *  ip ;
}
