#version 140

uniform mat4 model_view ;
uniform mat4 projection ;

in vec3 vertex ;
in float density ;
in mat4 frame ;
//varying vec3 pos ;

out vec3 normal ;
out vec3 pos ;
out float alpha ;

void main()
{
//  pos = vec3(gl_ModelViewMatrix * gl_Vertex);
//  normal = normalize( gl_NormalMatrix * gl_Vertex.xyz );

	normal = transpose( mat3(model_view) ) * vertex ;
	vec4 ip = frame * vec4( vertex, 1)  ;
	pos = ( model_view * ip ).xyz ;

	alpha = density ;

	gl_Position = projection * vec4( pos, 1) ;
}
