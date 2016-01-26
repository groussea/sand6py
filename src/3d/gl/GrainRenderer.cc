#include "GrainRenderer.hh"

#include "Texture.hh"
#include "ShapeRenderer.hh"

#include "visu/Offline.hh"

namespace d6 {


void GrainRenderer::draw_grains ( const Shader &shader, const float pixelSize,
								  const Eigen::Matrix4f &depthMVP, bool instanced ) const
{
	const int  divisor   = instanced ? 1 : 0 ;

	//Attributes
	gl::VertexAttribPointer vap_v( m_grainVertices, shader.attribute("vertex"), false, divisor ) ;
	gl::VertexAttribPointer vap_n( m_grainNormals, shader.attribute("normal"), false, divisor ) ;
	gl::VertexAttribPointer vap_a( m_grainVisibility, shader.attribute("visibility"), false, divisor ) ;

	// Uniforms
	glUniform1f( shader.uniform("grain_size"), m_offline.config().grainDiameter * m_grainSizeFactor ) ;
	glUniform1f( shader.uniform("pixel_size"), pixelSize ) ;

	glUniformMatrix4fv( shader.uniform("depth_mvp"), 1, GL_FALSE, depthMVP.data()) ;

	//vertices
	if( instanced )
	{
		gl::VertexPointer vp( m_shapeRenderer.squareVertices() ) ;
		glDrawArraysInstanced( GL_QUADS, 0, m_shapeRenderer.squareVertices().size(), m_grainVertices.size() );
	} else {

		if( pixelSize > 0 )
		{
			glEnable( GL_PROGRAM_POINT_SIZE ) ;
		} else {
			glPointSize( m_grainSizeFactor ) ;
		}

		gl::VertexPointer vp( m_grainVertices ) ;
		glDrawArrays( GL_POINTS, 0, m_grainVertices.size() );

		glDisable( GL_PROGRAM_POINT_SIZE ) ;
	}


}

void GrainRenderer::compute_shadow(
		const float pixelSize, const Eigen::Matrix4f &depthMVP )
{
	if( !m_depthShader.ok() )
		return ;

	UsingShader sh( m_depthShader ) ;

	const bool instanced = m_sampler.mode() == Sampler::Discs ;
	draw_grains( m_depthShader, pixelSize, depthMVP, instanced ) ;

}

void GrainRenderer::draw(const Texture& depthTexture, const Eigen::Vector3f &lightPosition,
		const float pixelSize, const Eigen::Matrix4f &depthMVP )
{
	if( !m_grainsShader.ok() )
		return ;

	UsingShader sh( m_grainsShader ) ;
	// Model-view
	sh.bindMVP() ;

	// Attributes
	const bool instanced = m_sampler.mode() == Sampler::Discs ;
	const int  divisor   = instanced ? 1 : 0 ;

	gl::VertexAttribPointer vap_s( m_grainNoise, m_grainsShader.attribute("noise"), false, divisor ) ;

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
}

} //d6
