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
		assert( valid() ) ;
		glDeleteBuffers( 1, &m_vbo ) ;
		m_vbo = INVALID_VBO ;
		m_size = 0 ;
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
		glVertexPointer(d, Traits< Scalar >::id , stride * sizeof( Scalar ),
						(GLvoid*) (offset * sizeof( Scalar )));
	}
	void set_normal_pointer( GLsizei stride = 0, GLuint offset = 0 ) const
	{
		bind() ;
		glNormalPointer( Traits< Scalar >::id , stride * sizeof( Scalar ),
						(GLvoid*)  (offset * sizeof( Scalar )));
	}
	void set_color_pointer( GLsizei stride = 0, GLuint offset = 0, unsigned d = 4 ) const
	{
		bind() ;
		glColorPointer(d, Traits< Scalar >::id , stride * sizeof( Scalar ),
					   (GLvoid*) (offset * sizeof( Scalar )));
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
		glEnableClientState( GL_VERTEX_ARRAY );
	}
	~VertexPointer()
	{
		glDisableClientState( GL_VERTEX_ARRAY );
	}
};
struct NormalPointer
{
	template< typename Scalar, unsigned Dim, int Type >
	NormalPointer( const VertexBuffer< Scalar, Dim, Type > &vb )
	{
		vb.set_normal_pointer();
		glEnableClientState( GL_NORMAL_ARRAY );
	}
	~NormalPointer()
	{
		glDisableClientState( GL_NORMAL_ARRAY );
	}
};
struct ColorPointer
{
	template< typename Scalar, unsigned Dim, int Type >
	ColorPointer( const VertexBuffer< Scalar, Dim, Type > &vb )
	{
		vb.set_color_pointer();
		glEnableClientState( GL_COLOR_ARRAY );
	}
	~ColorPointer()
	{
		glDisableClientState( GL_COLOR_ARRAY );
	}
};
struct VertexAttribPointer
{
	template< typename Scalar, unsigned Dim, int Type >
	VertexAttribPointer( const VertexBuffer< Scalar, Dim, Type > &vb, GLint attrib, 
		bool normalized = false )
		: m_attrib( attrib )
	{
		vb.set_vertex_attrib_pointer( attrib, normalized );
		glEnableVertexAttribArray( attrib );
	}
	~VertexAttribPointer()
	{
		glDisableVertexAttribArray( m_attrib );
	}
private:
	GLint m_attrib ;
};

typedef VertexBuffer< GLfloat, 4, GL_ARRAY_BUFFER > VertexBuffer4f ;
typedef VertexBuffer< GLfloat, 3, GL_ARRAY_BUFFER > VertexBuffer3f ;
typedef VertexBuffer< GLfloat,16, GL_ARRAY_BUFFER > VertexBuffer16f ;
typedef VertexBuffer< GLfloat, 1, GL_ARRAY_BUFFER > ArrayBufferf;

typedef VertexBuffer< GLdouble, 3, GL_ARRAY_BUFFER > VertexBuffer3d ;

typedef VertexBuffer< GLuint, 1, GL_ELEMENT_ARRAY_BUFFER > IndexBuffer ;
typedef VertexBuffer< GLuint, 1, GL_ARRAY_BUFFER > ArrayBufferui;

} //ns gl
} //ns d6

#endif
