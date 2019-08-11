/*
 * This file is part of Sand6, a C++ continuum-based granular simulator.
 *
 * Copyright 2016 Gilles Daviet <gilles.daviet@inria.fr> (Inria - Université Grenoble Alpes)
 *
 * Sand6 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * Sand6 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with Sand6.  If not, see <http://www.gnu.org/licenses/>.
*/

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

	void bindMVP( const float* modelView, const float* projection,
				  const char* model_view_name = "model_view",
				  const char* projection_name = "projection" ) ;
	void bindMVP( const char* model_view = "model_view",
				  const char* projection = "projection" ) ;

private:
	const Shader& m_shader ;
	GLuint previous ;
};

} //d6


#endif
