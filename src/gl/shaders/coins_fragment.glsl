#version 130

uniform mat4 model_view ;
uniform mat4 projection ;
uniform vec3 light_pos ;
uniform float grain_size ;
uniform sampler2D depth_texture;

in vec3 normal_screen ;
in vec3 vertex_screen ;
in vec4 shadow_coord  ;

in float vis ;
in float material ;

out vec4 color ;

void main (void)
{
    if( vis < 0 ) discard ;

    vec2 nproj = normalize( normal_screen.xy )  ;
    //Assume GL_POINT_SPRITE_COORD_ORIGIN is GL_UPPER_LEFT
    vec2 pos = vec2( gl_PointCoord.x - .5, .5 - gl_PointCoord.y ) ;
    vec2 nrot = vec2( -nproj[1], nproj[0] ) ;

    float h = (normal_screen[2] + 0.1)/1.1 ; //max( nproj[2] * nproj[2]) ;

    mat2 A  = 1./(h*h) * outerProduct(nproj.xy, nproj.xy) + outerProduct(nrot, nrot) ;

    if( dot( pos, A*pos) > 0.25 )
        discard ;


//    if( abs(dot(nproj.st, pos)) > .25 )
//        discard ;
//    if( gl_PointCoord.t > 0. )
//        discard ;

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

//        ambientMat += vec4( vec3(0.3,0.3, 0.), 1. );
//        diffuseMat += vec4( vec3(.6 , .6, 0.), 1. );
//        specMat    += vec4( alpha*alpha*vec3(.7 , .7, .7), 1. );

        ambientMat += mat_0*vec4( vec3(0.2,0.3, 0.1), 1. );
        diffuseMat += mat_0*vec4( vec3(.2 , .3, 0.1), 1. );

        ambientMat += mat_1*vec4( vec3(0.3,0.2, 0.1), 1. );
        diffuseMat += mat_1*vec4( vec3(.3 , .2, 0.1), 1. );

        ambientMat += mat_2*vec4( vec3(0.25,0.2, 0.15), 1. );
        diffuseMat += mat_2*vec4( vec3(.3, .2, 0.2), 1. );

        specMat    = vec4( 1., 1., 1., 1. );


    float specPow = 15.0;

    vec4 diffuse;
    vec4 spec;
    vec4 ambient;

    vec4 light = model_view * vec4( light_pos, 1 )  ;
    vec3 L = normalize( light.xyz - vertex_screen);
    vec3 E = normalize(- vertex_screen);
    vec3 R = normalize(reflect(-L, normal_screen ));

    ambient = ambientMat;
    diffuse = clamp( diffuseMat * max(abs(dot(normal_screen,L)), 0.0)  , 0.0, 1.0 ) ;
    spec = clamp ( specMat * pow(max(abs(dot(R,E)),0.0),0.3*specPow) , 0.0, 1.0 );

    color = ( 0.4 + 0.6*alpha) * ambient ;
    color += alpha * 0.75*diffuse ;
    color += alpha*alpha*0.5*spec;

    color.a = 0.5 ;
}
