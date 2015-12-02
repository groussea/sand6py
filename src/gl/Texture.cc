#include "Texture.hh"

#include "utils/Log.hh"

namespace d6 {



UsingFrameBuffer::UsingFrameBuffer( const FrameBuffer& fb ) {
	fb.bind();

	glGetIntegerv( GL_VIEWPORT, m_viewport );
	glViewport(0,0,fb.width(),fb.height()) ;
}

UsingFrameBuffer::~UsingFrameBuffer() {
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glViewport(0,0,m_viewport[2],m_viewport[3]) ;
}

UsingTexture::UsingTexture(const Texture &texture, GLuint unit)
	: m_texture(texture), m_unit(unit)
{
	glActiveTexture( unitID() ) ;
	texture.bind();
}

void UsingTexture::bindUniform(const GLint location) {
	glUniform1i( location, m_unit);
}

UsingTexture::~UsingTexture()
{
	glActiveTexture( unitID() ) ;
	glBindTexture( m_texture.target(), 0 ) ;
}

GLenum UsingTexture::unitID() const
{
	return GL_TEXTURE0 + m_unit ;
}

} //d6
