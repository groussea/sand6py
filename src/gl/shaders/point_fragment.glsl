

#version 330

layout(location = 0) out vec4 color ;

void main (void)
{
    vec4 ambientMat = vec4(0.3,0.2, 0.1, 1. );
    color = ambientMat ;
}

