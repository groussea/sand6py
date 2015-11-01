#version 110


varying vec3 normal;
varying vec3 pos ;
varying float alpha ;

void main (void)
{
	vec4 ambientMat = vec4(1.,0.,0., alpha );
	vec4 diffuseMat = vec4(0.,1.,0., alpha );
	vec4 specMat    = vec4( 0., 0., 1., alpha );
	float specPow = 15.0;

	vec4 diffuse;
	vec4 spec;
	vec4 ambient;

   vec3 L = normalize(gl_LightSource[0].position.xyz - pos);
   vec3 E = normalize(-pos);
   vec3 R = normalize(reflect(-L,normal));

	ambient = ambientMat;
	diffuse = clamp( diffuseMat * max(dot(normal,L), 0.0)  , 0.0, 1.0 ) ;
	spec = clamp ( specMat * pow(max(dot(R,E),0.0),0.3*specPow) , 0.0, 1.0 );

	gl_FragColor = ambient + diffuse + spec;
}
