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

#include "GrainRenderer.hh"

#include "Texture.hh"
#include "ShapeRenderer.hh"

#include "visu/Offline.hh"

namespace d6 {

void GrainRenderer::setup_vao( const Shader &shader, bool instanced)
{
	gl::ArrayObject::Using vao(m_grainArrays);

	//Attributes
	const int  divisor   = instanced ? 1 : 0 ;
	gl::VertexAttribPointer vap_v( m_grainVertices, shader.attribute("vertex"), false, divisor ) ;
	gl::VertexAttribPointer vap_n( m_grainNormals, shader.attribute("normal"), false, divisor ) ;
	gl::VertexAttribPointer vap_a( m_grainVisibility, shader.attribute("visibility"), false, divisor ) ;

	gl::VertexAttribPointer vap_s( m_grainNoise, m_grainsShader.attribute("noise"), false, divisor ) ;
}

void GrainRenderer::draw_grains ( const Shader &shader, const float pixelSize,
								  const Eigen::Matrix4f &depthMVP, bool instanced ) const
{

	gl::ArrayObject::Using vao(m_grainArrays);

	// Uniforms
	glUniform1f( shader.uniform("grain_size"), m_offline.config().grainDiameter * m_grainSizeFactor ) ;
	glUniform1f( shader.uniform("pixel_size"), pixelSize ) ;

	glUniformMatrix4fv( shader.uniform("depth_mvp"), 1, GL_FALSE, depthMVP.data()) ;

	//vertices
	if( instanced )
	{
		glDrawArraysInstanced( GL_QUADS, 0, m_shapeRenderer.squareVertices().size(), m_grainVertices.size() );
	} else {

		if( pixelSize > 0 )
		{
			glEnable( GL_PROGRAM_POINT_SIZE ) ;
		} else {
			glPointSize( m_grainSizeFactor ) ;
		}

		glDrawArrays( GL_POINTS, 0, m_grainVertices.size() );
		glDisable( GL_PROGRAM_POINT_SIZE ) ;
	}


}

void GrainRenderer::compute_shadow(
		const float pixelSize, const Eigen::Matrix4f &depthMVP ) const
{
	if( !m_depthShader.ok() )
		return ;

	UsingShader sh( m_depthShader ) ;

	const bool instanced = m_sampler.mode() == Sampler::Discs ;
	draw_grains( m_depthShader, pixelSize, depthMVP, instanced ) ;
}

void GrainRenderer::draw(const Texture& depthTexture, const Eigen::Vector3f &lightPosition,
		const Eigen::Matrix4f& modelView, const Eigen::Matrix4f& projection,
		const float pixelSize, const Eigen::Matrix4f &depthMVP ) const
{
	if( !m_grainsShader.ok() )
		return ;

	UsingShader sh( m_grainsShader ) ;
	// Model-view
	sh.bindMVP(modelView.data(), projection.data()) ;

	// Attributes
	const bool instanced = m_sampler.mode() == Sampler::Discs ;

	// Uniforms
	glUniform3fv( m_grainsShader.uniform("light_pos"), 1, lightPosition.data() ) ;

	UsingTexture tx( depthTexture ) ;
	tx.bindUniform( m_grainsShader.uniform("depth_texture") );

	draw_grains( m_grainsShader, pixelSize, depthMVP, instanced ) ;
}

void GrainRenderer::init()
{
	m_sampler.sampleParticles( m_nSamples ) ;

	m_grainsShader.add_attribute("vertex") ;
	m_grainsShader.add_attribute("normal") ;
	m_grainsShader.add_attribute("visibility") ;
	m_grainsShader.add_attribute("noise") ;

	m_grainsShader.add_uniform("model_view") ;
	m_grainsShader.add_uniform("projection") ;
	m_grainsShader.add_uniform("light_pos") ;
	m_grainsShader.add_uniform("depth_mvp");
	m_grainsShader.add_uniform("depth_texture");
	m_grainsShader.add_uniform("grain_size");
	m_grainsShader.add_uniform("pixel_size");

	m_depthShader.add_attribute("vertex") ;
	m_depthShader.add_attribute("normal") ;
	m_depthShader.add_attribute("visibility") ;
	m_depthShader.add_uniform("depth_mvp");
	m_depthShader.add_uniform("grain_size");
	m_depthShader.add_uniform("pixel_size");

	switch( m_sampler.mode() ){
	case Sampler::VelocityCut:
		m_grainsShader.load("grains_vertex","grains_vel_fragment") ;
		m_depthShader.load("grain_depth_vertex","grain_depth_fragment") ;
		break ;
	case Sampler::Discs:
		m_grainsShader.load("coins_vertex","coins_fragment") ;
		m_depthShader.load("coins_depth_vertex","coins_depth_fragment") ;
		break ;
	default:
		m_grainsShader.load("grains_vertex","grains_fragment") ;
		m_depthShader.load("grain_depth_vertex","grain_depth_fragment") ;
	}
}

void GrainRenderer::move()
{
	m_sampler.reassign();
	m_sampler.move( ) ;
}

void GrainRenderer::update_buffers()
{
	m_sampler.compute_absolute();

	m_grainVertices  .reset( m_sampler.count(), m_sampler.positions().data() , GL_DYNAMIC_DRAW )  ;
	m_grainNormals   .reset( m_sampler.count(), m_sampler.normals().data()   , GL_DYNAMIC_DRAW )  ;
	m_grainVisibility.reset( m_sampler.count(), m_sampler.visibility().data(), GL_DYNAMIC_DRAW )  ;
	m_grainNoise     .reset( m_sampler.count(), m_sampler.noise().data()     , GL_DYNAMIC_DRAW )  ;

	// Attributes
	const bool instanced = m_sampler.mode() == Sampler::Discs ;
	setup_vao(m_grainsShader, instanced);
}

} //d6
