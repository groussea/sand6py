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
//	draw() ; return ;
	{
		glPointSize( 3 );
		gl::VertexPointer vp( m_centers ) ;
		
		if( m_drawParticles ) {
			gl::ColorPointer  cp( m_colors ) ;
		} else {
			glColor3f(0,0,1) ;
		}
		
		glDrawArrays( GL_POINTS, 0, m_centers.size() );
	}

	for( const LevelSet::Ptr& ls: m_offline.levelSets() ) {
		drawObject( *ls );
	}

}

void GLViewer::draw()
{

	if( m_drawParticles ) 
	{
		m_glyphQuadIndices.bind();

		if( m_enableBending ) {
			glEnable (GL_BLEND);
			glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		}

		if( m_shader.ok() ) {

			UsingShader sh( m_shader ) ;

			// Model-view
			sh.bindMVP() ;

			//Vertices
			gl::VertexAttribPointer vap( m_glyph, m_shader.attributes["vertex"] ) ;

			// Densities
			gl::VertexAttribPointer  ap( m_alpha, m_shader.attributes["alpha"] ) ;
			glVertexAttribDivisor( m_shader.attributes["alpha"], 1 ) ;

			//Frames
			glEnableVertexAttribArray( m_shader.attributes["frame"]+0) ;
			glEnableVertexAttribArray( m_shader.attributes["frame"]+1) ;
			glEnableVertexAttribArray( m_shader.attributes["frame"]+2) ;
			glEnableVertexAttribArray( m_shader.attributes["frame"]+3) ;
			m_frames.set_vertex_attrib_pointer( m_shader.attributes["frame"]+0, false, 16, 4*0, 4 ) ;
			m_frames.set_vertex_attrib_pointer( m_shader.attributes["frame"]+1, false, 16, 4*1, 4 ) ;
			m_frames.set_vertex_attrib_pointer( m_shader.attributes["frame"]+2, false, 16, 4*2, 4 ) ;
			m_frames.set_vertex_attrib_pointer( m_shader.attributes["frame"]+3, false, 16, 4*3, 4 ) ;
			glVertexAttribDivisor( m_shader.attributes["frame"]+0, 1 ) ;
			glVertexAttribDivisor( m_shader.attributes["frame"]+1, 1 ) ;
			glVertexAttribDivisor( m_shader.attributes["frame"]+2, 1 ) ;
			glVertexAttribDivisor( m_shader.attributes["frame"]+3, 1 ) ;

			glDrawElementsInstanced( GL_QUADS, m_glyphQuadIndices.size(), GL_UNSIGNED_INT, 0, m_matrices.cols() );

			glDisableVertexAttribArray( m_shader.attributes["frame"]+3) ;
			glDisableVertexAttribArray( m_shader.attributes["frame"]+2) ;
			glDisableVertexAttribArray( m_shader.attributes["frame"]+1) ;
			glDisableVertexAttribArray( m_shader.attributes["frame"]+0) ;

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

		glDisable (GL_BLEND);
	}

	if( m_renderSamples )
	{
		if( m_grainsShader.ok() ) {

			UsingShader sh( m_grainsShader ) ;
			// Model-view
			sh.bindMVP() ;

			//Vertices
			gl::VertexAttribPointer vap_v( m_grainVertices, m_grainsShader.attribute("vertex") ) ;
			gl::VertexAttribPointer vap_n( m_grainNormals, m_grainsShader.attribute("normal") ) ;
		
			glPointSize( 2 ) ;
			glColor3f( 1,0,0 ) ;

			gl::VertexPointer vp( m_grainVertices ) ;
			glDrawArrays( GL_POINTS, 0, m_grainVertices.size() );
		
		}
	}

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
	Eigen::Matrix4f mat ;
	mat.setIdentity() ;
	mat.block<3,3>(0,0) = ls.rotation().matrix().cast < GLfloat >() * ls.scale()  ;
	mat.block<3,1>(0,3) = ls.origin().cast < GLfloat >() ;

	glPushMatrix();
	glMultMatrixf( mat.data() );

	glColor4f(1., .8, .8, 1);

	if( dynamic_cast<const SphereLevelSet*>(&ls) )
	{
		m_glyphQuadIndices.bind();
		gl::VertexPointer vp( m_glyph ) ;
//		gl::NormalPointer np( m_glyph ) ;
		glDrawElements( GL_QUADS, m_glyphQuadIndices.size(), GL_UNSIGNED_INT, 0 );
	}
	else if ( dynamic_cast<const PlaneLevelSet*>(&ls) )
	{
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


	m_shader.add_attribute("vertex") ;
	m_shader.add_attribute("frame") ;
	m_shader.add_attribute("alpha") ;

	m_shader.add_uniform("model_view") ;
	m_shader.add_uniform("projection") ;
	m_shader.load() ;

	if( m_renderSamples ) {
		m_sampler.sampleParticles( 10 ) ;

		m_grainsShader.add_attribute("vertex") ;
		m_grainsShader.add_attribute("normal") ;

		m_grainsShader.add_uniform("model_view") ;
		m_grainsShader.add_uniform("projection") ;
		m_grainsShader.load("grains_vertex","grains_fragment") ;
	}

	update_buffers();
}

void GLViewer::animate()
{
	if( ! next_frame() )
		stopAnimation();
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

	if( m_renderSamples ) {
		m_grainVertices.reset( m_sampler.count(), m_sampler.positions().data(), GL_DYNAMIC_DRAW )  ;
		m_grainNormals.reset( m_sampler.count(), m_sampler.normals().data(), GL_DYNAMIC_DRAW )  ;
	
	}
}

void GLViewer::postSelection(const QPoint& )
{
	Log::Info() << "Selected particle: " << selectedName() << std::endl ;
}


void GLViewer::set_frame(unsigned frame)
{
	if( frame == m_currentFrame+1 && m_renderSamples ) 
		m_sampler.updateOffsets( m_offline.frame_dt() * .5 ) ;

	if( m_offline.load_frame( frame ) ) {

		if( frame == m_currentFrame+1 && m_renderSamples ) 
			m_sampler.updateOffsets( m_offline.frame_dt() * .5 ) ;

		m_currentFrame = frame ;
	}

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
