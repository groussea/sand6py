#ifndef D6_SHADER_HH
#define D6_SHADER_HH

#include "opengl.hh"

#include <string>
#include <unordered_map>
#include <cassert>

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

	GLint uniform( const char* name ) const
	{
		const Bindings::const_iterator it = uniforms.find(name) ;
		assert( it != uniforms.end() ) ;
		return it->second ;
	}
	GLint attribute( const char* name ) const
	{
		const Bindings::const_iterator it = attributes.find(name) ;
		assert( it != attributes.end() ) ;
		return it->second ;
	}

};

struct UsingShader {
	UsingShader( const Shader& shader ) ;
	~UsingShader() ;

	void bindMVP( const char* model_view = "model_view",
				  const char* projection = "projection" ) ;

private:
	const Shader& m_shader ;
	GLuint previous ;
};

} //d6


#endif
