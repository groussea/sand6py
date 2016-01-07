#version 130

uniform mat4 model_view ;
uniform mat4 projection ;
uniform mat4 depth_mvp  ;
uniform float grain_size  ;

in vec3 vertex ;
in vec3 normal ;
in float visibility ;
in float noise ;

out vec2 coord ;
out vec3 normal_eye ;
out vec3 vertex_eye ;
out vec4 shadow_coord  ;
out float vis ;
out float material ;

void main()
{
    coord.y = 1 - 2*((gl_VertexID>>1)%2) ;
    coord.x = 1 - 2*(((gl_VertexID>>1)%2) ^ (gl_VertexID%2)) ;

    vec3 t1, t2 ;
    if( normal.x > 1.e-12 ) {
        t1 = normalize( vec3( -normal.y, normal.x, 0 ) ) ;
    } else {
        t1 = normalize( vec3( 0, -normal.z, normal.y ) ) ;
    }
    t2 = cross( normal, t1 );

    vec3 vertex_world = vertex + .5 * grain_size * ( coord.x*t1 + coord.y*t2 ) ;

    vertex_eye = ( model_view * vec4(vertex_world,1) ).xyz ;
    normal_eye  = normalize( model_view * vec4(normal,0) ).xyz ;

    gl_Position = projection * vec4(vertex_eye,1) ;
    shadow_coord = ( depth_mvp * vec4( vertex_world, 1 ) ) ;

    vis = visibility ;
    material = noise ;
}


