#version 330

uniform mat4 model_view ;
uniform mat4 projection ;

layout(location = 0) in vec3 vertex ;

void main()
{
	vec4 ip = model_view * vec4(vertex, 1)  ;

	gl_Position = projection * ip ;
}
