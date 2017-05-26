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

#include "Shader.hh"

#include "utils/Log.hh"
#include "utils/File.hh"

/*
 * Adapted from:
 * http://duriansoftware.com/joe/An-intro-to-modern-OpenGL.-Chapter-2.2:-Shaders.html
 */

namespace d6 {

#define show_info_log( object, glGet__iv, glGet__InfoLog ) \
{ \
	GLint log_length; \
	glGet__iv(object, GL_INFO_LOG_LENGTH, &log_length); \
	char * log = new char[log_length]; \
	glGet__InfoLog(object, log_length, NULL, log); \
	Log::Error() << log << std::endl ; \
	delete[] log; \
}

GLuint Shader::make_program(GLuint vertex_shader, GLuint fragment_shader)
{
	GLint program_ok;

	GLuint program = glCreateProgram();
	glAttachShader(program, vertex_shader);
	glAttachShader(program, fragment_shader);
	glLinkProgram(program);

	glGetProgramiv(program, GL_LINK_STATUS, &program_ok);
	 if (!program_ok) {
		 Log::Error() << "Failed to link shader program" << std::endl ;
		 show_info_log(program, glGetProgramiv, glGetProgramInfoLog);
		 glDeleteProgram(program);
		 return 0;
	 }
	 return program;
 }

Shader::Shader()
	: vertex_shader(0), fragment_shader(0), program(0)
{
}

bool Shader::load(const char *vertex, const char *fragment)
{
	vertex_shader = load( vertex, GL_VERTEX_SHADER ) ;
	fragment_shader = load( fragment, GL_FRAGMENT_SHADER ) ;

	if( vertex_shader && fragment_shader )
		program = make_program( vertex_shader, fragment_shader ) ;


	if( program ) {
		for( typename Bindings::iterator it = attributes.begin() ; it != attributes.end() ; ++ it ) {
			it->second = glGetAttribLocation( program, it->first.c_str() );
		}
		for( typename Bindings::iterator it = uniforms.begin() ; it != uniforms.end() ; ++ it ) {
			it->second = glGetUniformLocation( program, it->first.c_str() );
		}
	}

	return program != 0 ;
}

void Shader::destroy()
{
	if( program )
		glDeleteProgram( program );
	if( vertex_shader )
		glDeleteShader( vertex_shader );
	if( fragment_shader )
		glDeleteShader( fragment_shader );
}

Shader::~Shader()
{
	destroy() ;
}

GLuint Shader::load( const char* name, GLenum type )
{
	const std::string &dir =
			FileInfo(__FILE__).parentDirectory().filePath("shaders") ;

	const std::string filename = FileInfo( dir ).filePath( arg("%1.glsl", name) ) ;
	File file( filename, std::ios_base::in ) ;

	std::string source ;
	if( !file.get_contents( source ) )
	{
		Log::Error() << "Could not open shader file " << filename << std::endl ;
		return 0 ;
	}


	GLuint shader;
	shader = glCreateShader(type);

	const GLchar * source_c = source.c_str() ;
	GLint source_l = source.length() ;
	glShaderSource(shader, 1, &source_c, &source_l);

	glCompileShader(shader);

	GLint shader_ok;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &shader_ok);
	if (!shader_ok) {
		Log::Error() << arg( "Failed to compile %1:", filename) << std::endl ;
		show_info_log(shader, glGetShaderiv, glGetShaderInfoLog);
		glDeleteShader(shader);
		return 0;
	}

	return shader;
}

void Shader::use() const
{
	if( ok() ) {
		glUseProgram( program );
	}
}

UsingShader::UsingShader(const Shader &shader)
: m_shader( shader )
{
	glGetIntegerv( GL_CURRENT_PROGRAM, (GLint*) &previous );
	shader.use();
}

UsingShader::~UsingShader()
{
	glUseProgram( previous );
}

void UsingShader::bindMVP( const char* model_view_name, const char* projection_name )
{
	float model_view[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, model_view);
	float projection[16];
	glGetFloatv(GL_PROJECTION_MATRIX, projection);

	glUniformMatrix4fv(m_shader.uniform(model_view_name), 1, GL_FALSE, model_view );
	glUniformMatrix4fv(m_shader.uniform(projection_name), 1, GL_FALSE, projection );
}

} // d6
