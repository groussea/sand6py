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

#ifndef D6_TEXTURE_HH
#define D6_TEXTURE_HH

#include "opengl.hh"

#include <cassert>

namespace d6 {

struct Texture {

	Texture()
		:  m_id(INVALID_TEXTURE)
	{
	}

	~Texture()
	{
		destroy() ;
	}

	bool valid() const
	{
		return m_id != INVALID_TEXTURE ;
	}

	void gen()
	{
		assert( !valid() ) ;
		glGenTextures( 1, &m_id ) ;
	}

	void reset( GLenum target )
	{
		if( !valid() )
			gen() ;
		m_target = target ;
	}

	void bind( ) const
	{
		glBindTexture( m_target, m_id );
	}

	void destroy()
	{
		if( valid() ) {
			glDeleteTextures( 1, &m_id );
			m_id = INVALID_TEXTURE ;
		}
	}

	GLuint id() const {
		return m_id ;
	}

	GLenum target() const {
		return m_target ;
	}

private:
	static const GLuint INVALID_TEXTURE = -1 ;

	GLuint m_id ;
	GLenum m_target ;
};

struct FrameBuffer {

	FrameBuffer( )
		:  m_id(INVALID_FB), m_width(0), m_height(0)
	{}

	~FrameBuffer()
	{
		destroy() ;
	}

	bool valid() const
	{
		return m_id != INVALID_FB ;
	}

	void gen()
	{
		assert( !valid() ) ;
		glGenFramebuffers( 1, &m_id ) ;
	}

	void reset( GLuint width, GLuint height ) {
		if(!valid())
			gen() ;
		m_width = width ;
		m_height = height ;
	}

	void bind( GLenum target = GL_FRAMEBUFFER ) const
	{
		glBindFramebuffer( target, m_id );
	}

	void destroy()
	{
		if( valid() ) {
			glDeleteFramebuffers( 1, &m_id );
			m_id = INVALID_FB ;
		}
	}

	bool check_complete() {
		return glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE ;
	}

	GLuint width() const { return m_width ; }
	GLuint height() const { return m_height ; }

private:
	static const GLuint INVALID_FB = -1 ;

	GLuint m_id ;

	GLuint m_width ;
	GLuint m_height ;
};

struct UsingFrameBuffer
{

	UsingFrameBuffer( const FrameBuffer& fb )  ;
	~UsingFrameBuffer() ;

private:
	int m_viewport[4] ;
};

struct UsingTexture
{

	UsingTexture( const Texture& texture, GLuint unit = 0 )  ;

	void bindUniform( const GLint location ) ;

	~UsingTexture() ;

private:

	GLenum unitID() const ;

	const Texture & m_texture ;
	GLuint m_unit ;
};

} //d6

#endif
