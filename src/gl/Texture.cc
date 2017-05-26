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
