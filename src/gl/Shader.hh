#ifndef D6_SHADER_HH
#define D6_SHADER_HH

#include "opengl.hh"

#include <unordered_map>

namespace d6 {


class Shader {

public:
	Shader() ;
	~Shader() ;

	bool load( const char* vertex   = "vertex",
			   const char* fragment = "fragment"
			) ;
	void destroy() ;

	void use() const ;
	bool ok() const { return program != 0 ; }


	static GLuint load( const char* name, GLenum type ) ;
	static GLuint make_program(GLuint vertex_shader, GLuint fragment_shader) ;

	GLuint vertex_shader, fragment_shader, program;

	typedef std::unordered_map< std::string, GLint > Bindings ;
	Bindings attributes ;
	Bindings uniforms ;

	void add_attribute( const char* name ) 
	{ attributes[name] ; }
	void add_uniform( const char* name ) 
	{ uniforms[name] ; }

};

struct UsingShader {
	UsingShader( const Shader& shader ) ;
	~UsingShader() ;

private:
	GLuint previous ;
};

} //d6


#endif
