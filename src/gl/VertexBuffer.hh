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

#ifndef D6_GL_VERTEX_BUFFER_HH
#define D6_GL_VERTEX_BUFFER_HH

#include "opengl.hh"

#include <iostream>
#include <cassert>

namespace d6 {

namespace gl {

template < typename Scalar >
struct Traits
{
	static const unsigned id = GL_FLOAT ;
} ;
template <>
struct Traits< double >
{
	static const unsigned id = GL_DOUBLE ;
} ;

template< typename Scalar, unsigned Dim, int Type = GL_ARRAY_BUFFER  >
class VertexBuffer
{
public:

	VertexBuffer()
		:  m_vbo(INVALID_VBO), m_size( 0 )
	{
	}

	~VertexBuffer()
	{
		destroy() ;
	}

	bool valid() const
	{
		return m_vbo != INVALID_VBO ;
	}

	void gen()
	{
		assert( !valid() ) ;
		glGenBuffers( 1, &m_vbo ) ;
	}

	void destroy()
	{
		if( valid() ) {
			glDeleteBuffers( 1, &m_vbo ) ;
			m_vbo = INVALID_VBO ;
			m_size = 0 ;
		}
	}

	void bind() const
	{
		assert( valid() ) ;
//		std::cout << "Binding vbo " << m_vbo << std::endl ;
		glBindBuffer( Type, m_vbo ) ;
	}

	void set_vertex_pointer( GLsizei stride = 0, GLuint offset = 0, unsigned d = Dim  ) const
	{
		bind() ;
#ifndef GL_CORE
		glVertexPointer(d, Traits< Scalar >::id , stride * sizeof( Scalar ),
						(GLvoid*) (offset * sizeof( Scalar )));
#endif
	}
	void set_normal_pointer( GLsizei stride = 0, GLuint offset = 0 ) const
	{
		bind() ;
#ifndef GL_CORE
		glNormalPointer( Traits< Scalar >::id , stride * sizeof( Scalar ),
						(GLvoid*)  (offset * sizeof( Scalar )));
#endif
	}
	void set_color_pointer( GLsizei stride = 0, GLuint offset = 0, unsigned d = Dim ) const
	{
		bind() ;
#ifndef GL_CORE
		glColorPointer(d, Traits< Scalar >::id , stride * sizeof( Scalar ),
					   (GLvoid*) (offset * sizeof( Scalar )));
#endif
	}
	void set_vertex_attrib_pointer( GLint attribute, bool normalized = false, GLsizei stride = 0, GLuint offset = 0, unsigned d = Dim  ) const
	{
		bind() ;
		glVertexAttribPointer(attribute, d, Traits< Scalar >::id , normalized?GL_TRUE:GL_FALSE, stride * sizeof( Scalar ),
						(GLvoid*) (offset * sizeof( Scalar )));
	}

	void reset( GLuint size, const Scalar * data, GLenum usage = GL_STATIC_DRAW )
	{
		if( !valid() ) gen() ;
		bind() ;

		glBufferData(Type, size * Dim * sizeof(Scalar), data, usage );
		m_size = size ;

	}

	void update( GLuint size, const Scalar * data, GLuint offset = 0 )
	{
		bind() ;
		glBufferSubData(Type, offset * Dim * sizeof(Scalar), size * Dim * sizeof(Scalar), data );
	}

	template < typename WriteFunc >
	void map( const WriteFunc &func, GLenum mode = GL_WRITE_ONLY )
	{
		bind() ;
		Scalar * data = static_cast< Scalar * > (
					glMapBuffer( Type, mode ) ) ;
		func( data ) ;
		glUnmapBuffer( Type ) ;
	}

	GLuint size() const { return m_size ; }
	bool empty() const { return 0 == m_size ; }

private:
	static const GLuint INVALID_VBO = -1 ;

	GLuint m_vbo ;
	GLuint m_size ;

};

struct VertexPointer
{
	template< typename Scalar, unsigned Dim, int Type >
	VertexPointer( const VertexBuffer< Scalar, Dim, Type > &vb )
	{
		vb.set_vertex_pointer();
		#ifndef GL_CORE
		glEnableClientState( GL_VERTEX_ARRAY );
		#endif
	}
	~VertexPointer()
	{
		#ifndef GL_CORE
			glDisableClientState( GL_VERTEX_ARRAY );
		#endif
	}
};
struct NormalPointer
{
	template< typename Scalar, unsigned Dim, int Type >
	NormalPointer( const VertexBuffer< Scalar, Dim, Type > &vb )
	{
		if (vb.valid())
		{
			vb.set_normal_pointer();
#ifndef GL_CORE
			glEnableClientState(GL_NORMAL_ARRAY);
#endif
		}
	}
	~NormalPointer()
	{
		#ifndef GL_CORE
			glDisableClientState( GL_NORMAL_ARRAY );
		#endif
	}
};
struct ColorPointer
{
	template< typename Scalar, unsigned Dim, int Type >
	ColorPointer( const VertexBuffer< Scalar, Dim, Type > &vb )
	{
		if( vb.valid() ) {
			vb.set_color_pointer();
		#ifndef GL_CORE
			glEnableClientState( GL_COLOR_ARRAY );
#endif
		}
	}
	~ColorPointer()
	{
		#ifndef GL_CORE
		glDisableClientState( GL_COLOR_ARRAY );
		#endif
	}
};
struct VertexAttribPointer
{
	template< typename Scalar, unsigned Dim, int Type >
	VertexAttribPointer( const VertexBuffer< Scalar, Dim, Type > &vb, GLint attrib,
		bool normalized = false, int divisor = 0 )
		: m_attrib( attrib )
	{
		if( vb.valid() && vb.size() > 0 ) {
			vb.set_vertex_attrib_pointer( attrib, normalized );
			glEnableVertexAttribArray( attrib );
			glVertexAttribDivisor( m_attrib, divisor ) ;
		}
	}
	~VertexAttribPointer()
	{
		glVertexAttribDivisor( m_attrib,  0 ) ;
		glDisableVertexAttribArray( m_attrib );
	}
private:
	GLint m_attrib ;
};

template < unsigned Cols >
struct ArrayAttribPointer
{
	template< typename Scalar, unsigned Dim, int Type >
	ArrayAttribPointer( const VertexBuffer< Scalar, Dim, Type > &vb, GLint attrib,
		bool normalized = false, int divisor = 0 )
		: m_attrib( attrib )
	{
		for( unsigned i = 0 ; i < Cols ; ++i ) {
			glEnableVertexAttribArray( m_attrib+i );
			vb.set_vertex_attrib_pointer( m_attrib+i, normalized, Dim, i*Dim/Cols, Dim/Cols );
			glVertexAttribDivisor( m_attrib+i, divisor ) ;
		}
	}
	~ArrayAttribPointer()
	{
		for( unsigned i = 0 ; i < Cols ; ++i ) {
			glVertexAttribDivisor( m_attrib+i,  0 ) ;
			glDisableVertexAttribArray( m_attrib+i );
		}
	}
private:
	GLint m_attrib ;
};

struct VAO
{
	VAO()
		:  m_vao(INVALID_VAO)
	{
	}

	~VAO()
	{
		destroy() ;
	}

	bool valid() const
	{
		return m_vao != INVALID_VAO ;
	}

	void gen()
	{
		assert( !valid() ) ;
		glGenVertexArrays(1, &m_vao);
	}

	void destroy()
	{
		if( valid() ) {
			glDeleteVertexArrays( 1, &m_vao ) ;
			m_vao = INVALID_VAO ;
		}
	}

	void bind() const
	{
		glBindVertexArray(m_vao);
	}

	void unbind() const
	{
		glBindVertexArray(0);
	}

private:
	static const GLuint INVALID_VAO = -1 ;
	GLuint m_vao;
};

typedef VertexBuffer< GLfloat, 4, GL_ARRAY_BUFFER > VertexBuffer4f ;
typedef VertexBuffer< GLfloat, 3, GL_ARRAY_BUFFER > VertexBuffer3f ;
typedef VertexBuffer< GLfloat, 2, GL_ARRAY_BUFFER > VertexBuffer2f ;
typedef VertexBuffer< GLfloat,16, GL_ARRAY_BUFFER > VertexBuffer16f ;
typedef VertexBuffer< GLfloat, 1, GL_ARRAY_BUFFER > ArrayBufferf;

typedef VertexBuffer< GLdouble, 3, GL_ARRAY_BUFFER > VertexBuffer3d ;
typedef VertexBuffer< GLdouble, 2, GL_ARRAY_BUFFER > VertexBuffer2d ;

typedef VertexBuffer< GLuint, 1, GL_ELEMENT_ARRAY_BUFFER > IndexBuffer ;
typedef VertexBuffer< GLuint, 1, GL_ARRAY_BUFFER > ArrayBufferui;

} //ns gl
} //ns d6

#endif
