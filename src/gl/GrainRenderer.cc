#include "GrainRenderer.hh"

#include "Texture.hh"

#include "visu/Offline.hh"

namespace d6 {

void GrainRenderer::compute_shadow(
		const float pixelSize, const Eigen::Matrix4f &depthMVP )
{
	if( !m_grainsShader.ok() )
		return ;

	UsingShader sh( m_depthShader ) ;

	if( pixelSize > 0 )
	{
		glEnable( GL_PROGRAM_POINT_SIZE ) ;
		glUniform1f( m_depthShader.uniform("grain_size"), m_offline.config().grainDiameter * m_grainSizeFactor ) ;
		glUniform1f( m_depthShader.uniform("pixel_size"), pixelSize ) ;
	} else {
		glPointSize( 2 * m_grainSizeFactor ) ;
	}


	gl::VertexPointer vp( m_grainVertices ) ;

	glUniformMatrix4fv( m_depthShader.uniform("depth_mvp"), 1, GL_FALSE, depthMVP.data()) ;

	gl::VertexAttribPointer vap_v( m_grainVertices, m_depthShader.attribute("vertex") ) ;
	gl::VertexAttribPointer vap_a( m_grainVisibility, m_grainsShader.attribute("visibility") ) ;

	glDrawArrays( GL_POINTS, 0, m_grainVertices.size() );

	glDisable( GL_PROGRAM_POINT_SIZE ) ;

}

void GrainRenderer::draw(const Texture& depthTexture, const Eigen::Vector3f &lightPosition,
		const float pixelSize, const Eigen::Matrix4f &depthMVP )
{
	if( !m_grainsShader.ok() )
		return ;

	UsingShader sh( m_grainsShader ) ;
	// Model-view
	sh.bindMVP() ;

	if( m_sampler.mode() == Sampler::Discs )
	{
		glEnable( GL_POINT_SPRITE ) ;
	}
	if( pixelSize > 0 )
	{
		glEnable( GL_PROGRAM_POINT_SIZE ) ;
		glUniform1f( m_grainsShader.uniform("grain_size"), m_offline.config().grainDiameter * m_grainSizeFactor ) ;
		glUniform1f( m_grainsShader.uniform("pixel_size"), pixelSize ) ;
	} else {
		glPointSize( m_grainSizeFactor ) ;
	}

	//Vertices

	gl::VertexPointer vp( m_grainVertices ) ;
	gl::VertexAttribPointer vap_v( m_grainVertices, m_grainsShader.attribute("vertex") ) ;
	gl::VertexAttribPointer vap_n( m_grainNormals, m_grainsShader.attribute("normal") ) ;

	gl::VertexAttribPointer vap_a( m_grainVisibility, m_grainsShader.attribute("visibility") ) ;
	gl::VertexAttribPointer vap_s( m_grainNoise, m_grainsShader.attribute("noise") ) ;

	glUniform3fv( m_grainsShader.uniform("light_pos"), 1, lightPosition.data() ) ;

	UsingTexture tx( depthTexture ) ;
	tx.bindUniform( m_grainsShader.uniform("depth_texture") );

	glUniformMatrix4fv( m_grainsShader.uniform("depth_mvp"), 1, GL_FALSE, depthMVP.data()) ;

	glDrawArrays( GL_POINTS, 0, m_grainVertices.size() );

	glDisable( GL_PROGRAM_POINT_SIZE ) ;
	glDisable( GL_POINT_SPRITE ) ;

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
	if( m_sampler.mode() == Sampler::VelocityCut ) {
		m_grainsShader.load("grains_vertex","grains_vel_fragment") ;
	} else if( m_sampler.mode() == Sampler::Discs ) {
		m_grainsShader.load("grains_vertex","coins_fragment") ;
	} else {
		m_grainsShader.load("grains_vertex","grains_fragment") ;
	}

	m_depthShader.add_attribute("vertex") ;
	m_depthShader.add_uniform("depth_mvp");
	m_depthShader.add_uniform("grain_size");
	m_depthShader.add_uniform("pixel_size");
	m_depthShader.load("grain_depth_vertex","grain_depth_fragment") ;
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
