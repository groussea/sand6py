#version 130

uniform mat4 model_view ;
uniform vec3 light_pos ;
uniform sampler2D depth_texture;

in vec3 normal_eye ;
in vec3 vertex_eye ;
in vec4 shadow_coord  ;

in float vis ;
in float material ;

out vec4 color ;

void main (void)
{
    if( vis < 0 ) discard ;

    float alpha = 1 ; //vis ;

    vec2 tex_coords = 0.5 * ( shadow_coord.xy / shadow_coord.w ) + vec2(0.5,0.5) ;

    float zz = shadow_coord.z / shadow_coord.w  ;
    float zs = texture( depth_texture, tex_coords ).r ;

    if ( zs  <  zz ){
         alpha = 1. - (zz-zs)*10 ; //shadow_coord.z/20;
     }

    alpha = pow( clamp(alpha,0,1), 3) ;
    //alpha = vis ;
    color = alpha * vec4( vec3(0.3,0.3, 0.), 1. );
    //return ;

    vec4 ambientMat  = vec4(0,0,0,0) ;
    vec4 diffuseMat  = vec4(0,0,0,0) ;
    vec4 specMat     = vec4(0,0,0,0) ;

    float mat_0 = clamp( 1 - 2*(material-0.0), 0.0, 1.0 ) ;
    float mat_1 = clamp( 1 - 2*(material-0.5), 0.0, 1.0 ) ;
    float mat_2 = clamp( 1 - 2*(material-1.0), 0.0, 1.0 ) ;

    ambientMat += mat_0*vec4( vec3(0.2,0.3, 0.1), 1. );
    diffuseMat += mat_0*vec4( vec3(.2 , .3, 0.1), 1. );

    ambientMat += mat_1*vec4( vec3(0.3,0.2, 0.1), 1. );
    diffuseMat += mat_1*vec4( vec3(.3 , .2, 0.1), 1. );

    ambientMat += mat_2*vec4( vec3(0.25,0.2, 0.15), 1. );
    diffuseMat += mat_2*vec4( vec3(.3, .2, 0.2), 1. );

    vec4 diffuse;
    vec4 ambient;

    vec4 light = model_view * vec4( light_pos, 1 )  ;
    vec3 L = normalize( light.xyz - vertex_eye);

    diffuse = clamp( diffuseMat * max(dot(normal_eye,L), 0.0)  , 0.0, 1.0 ) ;

    color = 0.5 * ambientMat ;
    color += alpha * diffuse;

    color.a = 0.5 ;
}
