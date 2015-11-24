#version 130

uniform mat4 model_view ;
uniform vec3 light_pos ;
uniform sampler2D depth_texture;

in vec3 normal_screen ;
in vec3 vertex_screen ;
in vec4 shadow_coord  ;

in float vis ;
in float material ;

out vec4 color ;

void main (void)
{
    float alpha = 1 ; //vis ;

    vec2 tex_coords = 0.5 * ( shadow_coord.xy / shadow_coord.w ) + vec2(0.5,0.5) ;

    float zz = shadow_coord.z / shadow_coord.w  ;
    float zs = texture( depth_texture, tex_coords ).r ;

    if ( zs  <  zz ){
         alpha = 1. - shadow_coord.w*(zz-zs) ; //shadow_coord.z/20;
     }

    alpha = 0.6 * pow( clamp(alpha,0,1), 3) ;

    vec4 ambientMat  = vec4(0,0,0,0) ;
    vec4 diffuseMat  = vec4(0,0,0,0) ;
    vec4 specMat     = vec4(0,0,0,0) ;

    float mat_0 = clamp( 1 - 2*(material-0.0), 0.0, 1.0 ) ;
    float mat_1 = clamp( 1 - 2*(material-0.5), 0.0, 1.0 ) ;
    float mat_2 = clamp( 1 - 2*(material-1.0), 0.0, 1.0 ) ;

        ambientMat += mat_0*vec4( vec3(0.3,0.3, 0.), 1. );
        diffuseMat += mat_0*vec4( alpha*vec3(.8 , .8, 0.), 1. );
        specMat    += mat_0*vec4( alpha*alpha*vec3(.7 , .7, .7), 1. );

        ambientMat += mat_1*vec4( vec3(0.3,0.2, 0.), 1. );
        diffuseMat += mat_1*vec4( alpha*vec3(.8 , .6, 0.), 1. );
        specMat    += mat_1*vec4( alpha*alpha*vec3(.7 , .7, .7), 1. );

        ambientMat += mat_2*vec4( vec3(0.,0., 0.), 1. );
        diffuseMat += mat_2*vec4( alpha*vec3(.5, .5, 0.5), 1. );
        specMat    += mat_2*vec4( alpha*alpha*vec3(.7 , .7, .7), 1. );


    float specPow = 15.0;

    vec4 diffuse;
    vec4 spec;
    vec4 ambient;

    vec4 light = model_view * vec4( light_pos, 1 )  ;
    vec3 L = normalize( light.xyz - vertex_screen);
    vec3 E = normalize(- vertex_screen);
    vec3 R = normalize(reflect(-L, normal_screen ));

    ambient = ambientMat;
    diffuse = clamp( diffuseMat * max(dot(normal_screen,L), 0.0)  , 0.0, 1.0 ) ;
    spec = clamp ( specMat * pow(max(dot(R,E),0.0),0.3*specPow) , 0.0, 1.0 );

    color = ambient + diffuse + spec;
    color.a = 0.5 ;
}
