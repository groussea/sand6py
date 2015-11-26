#include "GLViewer.hh"

#include "visu/Offline.hh"

#include "geo/Particles.hh"
#include "geo/Tensor.hh"
#include "geo/Grid.hh"
#include "geo/LevelSet.impl.hh"

#include "utils/Log.hh"
#include "utils/File.hh"

#include <Eigen/Eigenvalues>

#include <GL/gl.h>

#include <QKeyEvent>
#include <QApplication>

#include <iomanip>

#define fb_width  1024
#define fb_height 1024

namespace d6 {

static void genSphere( const unsigned parallels, const unsigned meridians,
					   Eigen::Matrix3Xf& ballVerts,
					   std::vector< GLuint >& quadIndices )
{
	ballVerts.resize( 3, parallels * meridians ) ;

	quadIndices.clear();
	quadIndices.reserve( 4 * (parallels - 1) * meridians ) ;

	const double dphi = M_PI / (parallels - 1.);
	const double dpsi = 2 * M_PI / meridians ;

	for(unsigned i = 0; i < parallels; ++i)
	{
		const double z = std::cos( i * dphi ) ;
		const double r = std::sin( i * dphi ) ;

		for(unsigned j = 0; j < meridians; j++)
		{
			const Vec p ( r * std::cos( j * dpsi ), r * std::sin( j * dpsi ), z )  ;
			ballVerts.col( i*meridians + j ) = p.cast< GLfloat >() ;

			if( i > 0 ) {
				quadIndices.push_back( (i-1)*meridians + ( ( j+1 ) % meridians ) ) ;
				quadIndices.push_back( (i-1)*meridians + j ) ;
				quadIndices.push_back( i*meridians + j ) ;
				quadIndices.push_back( i*meridians +  ( ( j+1 ) % meridians ) ) ;
			}
		}
	}

}

void GLViewer::fastDraw()
{
	if( !m_fastDraw ) {
		draw() ;
		return ;
	}

	{
		glPointSize( 3 );
		gl::VertexPointer vp( m_centers ) ;

		if( m_drawParticles ) {
			gl::ColorPointer  cp( m_colors ) ;
			glDrawArrays( GL_POINTS, 0, m_centers.size() );
		} else {
			glColor3f(0,0,1) ;
			glDrawArrays( GL_POINTS, 0, m_centers.size() );
		}

	}

	for( const LevelSet::Ptr& ls: m_offline.levelSets() ) {
		drawObject( *ls );
	}

}

void GLViewer::draw()
{

	if( m_enableBending ) {
		glEnable (GL_BLEND);
		glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}


	if( m_drawParticles )
	{
		m_glyphQuadIndices.bind();

		if( m_shader.ok() ) {

			UsingShader sh( m_shader ) ;

			// Model-view
			sh.bindMVP() ;

			//Vertices
			gl::VertexAttribPointer vap( m_glyph, m_shader.attributes["vertex"] ) ;

			// Densities
			gl::VertexAttribPointer  ap( m_alpha, m_shader.attributes["alpha"], false, 1 ) ;

			//Frames
			gl::ArrayAttribPointer<4>  fp( m_frames, m_shader.attributes["frame"], false, 1 ) ;

			glDrawElementsInstanced( GL_QUADS, m_glyphQuadIndices.size(), GL_UNSIGNED_INT, 0, m_matrices.cols() );

		} else {

			gl::VertexPointer vp( m_glyph ) ;
			gl::NormalPointer np( m_glyph ) ;

			for( int i = 0 ; i < m_matrices.cols() ; ++i ){
				glPushMatrix();
				glMultMatrixf( m_matrices.col(i).data() );

				glColor4f(1., 0, 0, m_densities[i]);
				glDrawElements( GL_QUADS, m_glyphQuadIndices.size(), GL_UNSIGNED_INT, 0 );

				glPopMatrix();
			}

		}

	}

	if( renderSamples() )
	{
		if( m_grainsShader.ok() ) {

			Eigen::Matrix4f depthMVP ;

//			qglviewer::Camera& cam = *camera() ;
			qglviewer::Camera  cam = *camera() ;
			Eigen::Vector3f light_pos = lightPosition() ;
			qglviewer::Vec lp( light_pos[0], light_pos[1], light_pos[2] )  ;

			cam.setPosition( lp );
			cam.lookAt( sceneCenter() );
			cam.setFieldOfView( std::atan2( m_offline.mesh().box().norm(), lightPosition().norm()/2 ) ) ;
//			cam.setFieldOfView( 2./3*M_PI );
//			cam.setSceneRadius(  );
			cam.computeModelViewMatrix();
			cam.computeProjectionMatrix();

			Eigen::Matrix4d depthMVP_d ;
			cam.getModelViewProjectionMatrix(depthMVP_d.data());
			depthMVP = depthMVP_d.cast< float >() ;

			gl::VertexPointer vp( m_grainVertices ) ;

			if(1){

				glBindFramebuffer(GL_FRAMEBUFFER, m_depthBuffer );
				glBindTexture(GL_TEXTURE_2D, m_depthTexture);

				glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, m_depthTexture, 0);
				glDrawBuffer(GL_NONE); // No color buffer is drawn to.

				if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
					Log::Error() << "Frame buffer incomplete" << std::endl ;

				int viewport[4] ;
				glGetIntegerv( GL_VIEWPORT, viewport );
				glViewport(0,0,fb_width,fb_height) ;

				UsingShader sh( m_depthShader ) ;

				glUniformMatrix4fv( m_depthShader.uniform("depth_mvp"), 1, GL_FALSE, depthMVP.data()) ;

				gl::VertexAttribPointer vap_v( m_grainVertices, m_depthShader.attribute("vertex") ) ;

				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

				glPointSize( 2 ) ;
				glDrawArrays( GL_POINTS, 0, m_grainVertices.size() );

				glBindFramebuffer(GL_FRAMEBUFFER, 0);
				glViewport(0,0,viewport[2],viewport[3]) ;

//				*camera() = old_cam ;
			}

			if(0){
				UsingShader sh( m_testShader ) ;

				glActiveTexture(GL_TEXTURE0) ;
				glBindTexture( GL_TEXTURE_2D, m_depthTexture ) ;
				glUniform1i( m_testShader.uniform("in_texture"), 0);

				gl::VertexAttribPointer vap_v( m_square, m_testShader.attribute("vertex") ) ;
				gl::VertexPointer vp( m_square ) ;
				glDrawArrays( GL_QUADS, 0, m_square.size() ) ;

				glBindTexture( GL_TEXTURE_2D, 0 ) ;
			}

			if(1){

				UsingShader sh( m_grainsShader ) ;
				// Model-view
				sh.bindMVP() ;

				//Vertices
				gl::VertexAttribPointer vap_v( m_grainVertices, m_grainsShader.attribute("vertex") ) ;
				gl::VertexAttribPointer vap_n( m_grainNormals, m_grainsShader.attribute("normal") ) ;

				gl::VertexAttribPointer vap_a( m_grainVisibility, m_grainsShader.attribute("visibility") ) ;
				gl::VertexAttribPointer vap_s( m_grainNoise, m_grainsShader.attribute("noise") ) ;

				glUniform3fv( m_grainsShader.uniform("light_pos"), 1, lightPosition().data() ) ;

				glActiveTexture(GL_TEXTURE0) ;
				glBindTexture( GL_TEXTURE_2D, m_depthTexture ) ;
				glUniform1i( m_grainsShader.uniform("depth_texture"), 0);

				glUniformMatrix4fv( m_grainsShader.uniform("depth_mvp"), 1, GL_FALSE, depthMVP.data()) ;

				glPointSize( 1 ) ;
				glDrawArrays( GL_POINTS, 0, m_grainVertices.size() );

				glBindTexture( GL_TEXTURE_2D, 0 ) ;

			}

		}
	}

	glDisable (GL_BLEND);

	for( const LevelSet::Ptr& ls: m_offline.levelSets() ) {
		drawObject( *ls );
	}

	if( m_snapshotting )
		snap() ;

//	Log::Debug() << "Current fps " << currentFPS() << std::endl ;

}

void GLViewer::drawWithNames()
{

	if( m_drawParticles )
	{

		m_glyphQuadIndices.bind();

		gl::VertexPointer vp( m_glyph ) ;

		for( int i = 0 ; i < m_matrices.cols() ; ++i ){
			glPushMatrix();
			glMultMatrixf( m_matrices.col(i).data() );
			glPushName(i) ;

			glDrawElements( GL_QUADS, m_glyphQuadIndices.size(), GL_UNSIGNED_INT, 0 );

			glPopName() ;
			glPopMatrix();
		}

	}
}

void GLViewer::drawObject(const LevelSet &ls)
{
	const Eigen::Matrix3f rotation = ls.rotation().matrix().cast < GLfloat >() ;
	const Eigen::Vector3f translation = ls.origin().cast < GLfloat >() ;

	if( dynamic_cast<const SphereLevelSet*>(&ls) )
	{

		UsingShader sh( m_ballShader ) ;
		// Model-view
		sh.bindMVP() ;

		//Vertices
		gl::VertexAttribPointer vap_v( m_square, m_ballShader.attribute("vertex") ) ;

		glUniform1f( m_ballShader.uniform("radius"), ls.scale() ) ;
		glUniformMatrix3fv( m_ballShader.uniform("rotation"), 1, GL_FALSE, rotation.data() ) ;
		glUniform3fv( m_ballShader.uniform("center"), 1, translation.data() ) ;

		glUniform3fv( m_ballShader.uniform("light_pos"), 1, lightPosition().data() ) ;

		gl::VertexPointer vp( m_square ) ;
		glDrawArrays( GL_QUADS, 0, m_square.size() ) ;

	} else {

		Eigen::Matrix4f mat ;
		mat.setIdentity() ;
		mat.block<3,3>(0,0) = rotation * ls.scale()  ;
		mat.block<3,1>(0,3) = translation ;

		glColor4f(1., .8, .8, 1);


		glPushMatrix();
		glMultMatrixf( mat.data() );

		if ( dynamic_cast<const PlaneLevelSet*>(&ls) ) {
			const Vec& box = m_offline.mesh().box() ;
			glBegin( GL_QUADS );
			glNormal3f( 0.f, 0.f, 1.f );
			glVertex3d( -box[0], -box[1], 0 );
			glNormal3f( 0.f, 0.f, 1.f );
			glVertex3d( -box[0],  box[1], 0 );
			glNormal3f( 0.f, 0.f, 1.f );
			glVertex3d(  box[0],  box[1], 0 );
			glNormal3f( 0.f, 0.f, 1.f );
			glVertex3d(  box[0], -box[1], 0 );
			glEnd( ) ;
		}
		glPopMatrix();
	}


}

void GLViewer::init()
{
	// Restore previous viewer state.
	restoreStateFromFile();

	// Camera
	const Vec& box = m_offline.mesh().box() ;
	const qglviewer::Vec qgl_box( box[0], box[1], box[2] ) ;
	const qglviewer::Vec qgl_ori(0,0,0) ;
	setSceneBoundingBox( qgl_ori, qgl_box) ;
	camera()->setZClippingCoefficient( box.maxCoeff() );

	// Gen glyph vertices
	Eigen::Matrix3Xf sphereVertices ;
	std::vector< GLuint > quadIndices ;
	genSphere( 5, 8, sphereVertices, quadIndices );
	m_glyph.reset( sphereVertices.cols(), sphereVertices.data(), GL_STATIC_DRAW );
	m_glyphQuadIndices.reset( quadIndices.size(), quadIndices.data() );

	Eigen::Matrix<float, 3, 4> vtx ;
	vtx  <<  -1, -1, 1,  1,
		 -1,  1, 1, -1,
		 0, 0 ,0, 0 ;
	m_square.reset( 4, vtx.data() );


	m_shader.add_attribute("vertex") ;
	m_shader.add_attribute("frame") ;
	m_shader.add_attribute("alpha") ;

	m_shader.add_uniform("model_view") ;
	m_shader.add_uniform("projection") ;
	m_shader.load() ;

	if( renderSamples() ) {
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
		m_grainsShader.load("grains_vertex","grains_fragment") ;

		glGenFramebuffers(1, &m_depthBuffer);

		glGenTextures(1, &m_depthTexture);
		glBindTexture(GL_TEXTURE_2D, m_depthTexture);

		//float data[ fb_width * fb_height ]	;
		//for( unsigned j = 0 ; j < fb_height ; ++j )
		//	for( unsigned i = 0 ; i < fb_width ; ++i ) {
		//		daaa[ j*fb_width + i ] = (1.*i)/fb_width ;
		//	}
		float* data = 0 ;

		glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, fb_width, fb_height, 0,GL_DEPTH_COMPONENT, GL_FLOAT, data);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

		m_depthShader.add_attribute("vertex") ;
		m_depthShader.add_uniform("depth_mvp");
		m_depthShader.load("depth_vertex","depth_fragment") ;
	}

	m_ballShader.add_uniform("model_view") ;
	m_ballShader.add_uniform("projection") ;
	m_ballShader.add_uniform("radius") ;
	m_ballShader.add_uniform("rotation") ;
	m_ballShader.add_uniform("center") ;
	m_ballShader.add_uniform("light_pos") ;
	m_ballShader.add_attribute("vertex") ;
	m_ballShader.load("ball_vertex","ball_fragment") ;

	m_testShader.add_attribute("vertex") ;
	m_testShader.add_uniform("in_texture");
	m_testShader.load("textest_vertex","textest_fragment") ;

	update_buffers();
}

void GLViewer::animate()
{
	if( ! next_frame() )
		stopAnimation();
}

Eigen::Vector3f GLViewer::lightPosition() const
{
	return Eigen::Vector3f ( 0, 0, m_offline.mesh().box()[2] * 2 ) ;
}

void GLViewer::update_buffers()
{
	const Particles &p = m_offline.particles() ;
	m_centers.reset( p.count(), p.centers().data(), GL_STATIC_DRAW )  ;

	if( m_drawParticles )
	{

		m_matrices.resize( 16, p.count() );
		m_densities.resize( p.count() );

		// Compute movel-view matrix from tensor
		Eigen::Matrix4f mat ;
		Mat tensor ;

#pragma omp parallel for private(mat, tensor)
		for( size_t i = 0 ; i < p.count() ; ++i ) {
			mat.setIdentity()  ;

			tensor_view( p.frames().col( i ) ).get( tensor ) ;
			Eigen::SelfAdjointEigenSolver<Mat> es( tensor );

			const Vec ev = es.eigenvalues().array().max(0).sqrt() ;
			const Scalar vol = 8 * ev.prod() ;

			mat.block<3,3>(0,0) = ( es.eigenvectors() * ev.asDiagonal() ).cast< GLfloat >()  ;
			mat.block<3,1>(0,3) = p.centers().col(i).cast < GLfloat >() ;

			m_matrices.col(i) = Eigen::Matrix< GLfloat, 16, 1 >::Map( mat.data(), mat.size() ) ;
			m_densities[i] = p.volumes()[i] / vol ;

		}
		m_frames.reset( p.count(), m_matrices.data(), GL_STATIC_DRAW )  ;
		m_alpha.reset ( p.count(), m_densities.data(), GL_STATIC_DRAW )  ;

		// Colors
		Eigen::Matrix4Xf colors( 4, p.count() ) ;
		colors.topRows(2).setZero() ;
		colors.row(2) = m_densities.cast< float >() ;
		colors.row(3).setOnes() ;
		m_colors.reset( p.count(), colors.data(), GL_STATIC_DRAW )  ;
	}

	if( renderSamples() ) {
		m_grainVertices  .reset( m_sampler.count(), m_sampler.positions().data() , GL_DYNAMIC_DRAW )  ;
		m_grainNormals   .reset( m_sampler.count(), m_sampler.normals().data()   , GL_DYNAMIC_DRAW )  ;
		m_grainVisibility.reset( m_sampler.count(), m_sampler.visibility().data(), GL_DYNAMIC_DRAW )  ;
		m_grainNoise     .reset( m_sampler.count(), m_sampler.noise().data()     , GL_DYNAMIC_DRAW )  ;
	}
}

void GLViewer::postSelection(const QPoint& )
{
	Log::Info() << "Selected particle: " << selectedName() << std::endl ;
}


void GLViewer::set_frame(unsigned frame)
{

	if( frame == m_currentFrame+1 && renderSamples() ) {
		m_sampler.reassign();
		m_sampler.move( m_offline.frame_dt() ) ;
	}

	if( m_offline.load_frame( frame ) ) {

		m_currentFrame = frame ;
	}

	m_sampler.compute_absolute() ;

	update_buffers();
}

void GLViewer::keyPressEvent(QKeyEvent *e)
{
	switch (e->key())
	{
	case Qt::Key_I :
	case Qt::Key_Period :
		next_frame();
		break ;
	case Qt::Key_P :
		prev_frame();
		break ;
	case Qt::Key_Home :
		set_frame(0);
		break ;
	case Qt::Key_D:
		m_drawParticles = !m_drawParticles ;
		update_buffers() ;
		break ;
	case Qt::Key_F:
		m_fastDraw = !m_fastDraw ;
		break ;
	case Qt::Key_B:
		m_enableBending = !m_enableBending ;
		break ;
	case Qt::Key_R:
		m_snapshotting = !m_snapshotting ;
		break ;
	case Qt::Key_Q :
		QApplication::exit( 0 ) ;
	default:
		QGLViewer::keyPressEvent(e);
	}
	updateGL() ;
}

void GLViewer::snap() {
	if( m_lastSnapped != m_currentFrame ) {

		//Snaps
		FileInfo snap_dir( FileInfo( m_offline.base_dir() ).filePath("snaps") ) ;
		FileInfo snap_file( snap_dir.filePath("qgl-%1.png") ) ;
		snap_file.makePath() ;

		std::stringstream num ;
		num << std::setfill('0') << std::setw(4) <<  m_currentFrame ;
		const std::string fileName = arg(snap_file.path(),  num.str() ) ;

		QImage snapshot = grabFrameBuffer( true );
		bool ok = snapshot.save( QString( fileName.c_str() ), "PNG" );
		if( ok ) {
			Log::Verbose() << "Saved snap of frame " << m_currentFrame << " on " << fileName << std::endl ;
		} else {
			Log::Error() << "Failed snapping frame " << m_currentFrame << " on " << fileName << std::endl ;
		}
		m_lastSnapped = m_currentFrame ;
	}
}

} //d6
