#version 130

uniform mat4 model_view ;
uniform mat4 projection ;
uniform mat4 depth_mvp  ;
uniform float grain_size  ;
uniform float pixel_size  ;

in vec3 vertex ;
in vec3 normal ;
in float visibility ;
in float noise ;

out vec3 normal_screen ;
out vec3 vertex_screen ;
out vec4 shadow_coord  ;
out float vis ;
out float material ;

void main()
{
    vec4 ip = model_view * vec4( vertex, 1)  ;
    vertex_screen  = ip.xyz ;
    normal_screen  = normalize( mat3(model_view) * normal ) ;

    vis = visibility ;

    material = noise ;

    gl_Position = projection *  ip ;
    shadow_coord = ( depth_mvp * vec4( vertex, 1 ) ) ;

    gl_PointSize = max( 1, pixel_size * grain_size / gl_Position.w );
}
