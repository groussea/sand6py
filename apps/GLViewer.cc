#include "GLViewer.hh"

#include "visu/Offline.hh"
#include "geo/Particles.hh"

#include <GL/gl.h>

#include <QKeyEvent>
#include <QApplication>

namespace d6 {

void GLViewer::fastDraw()
{
	{
		glPointSize( 3 );
		gl::VertexPointer vp( m_centers ) ;
		gl::ColorPointer  cp( m_colors ) ;
		glDrawArrays( GL_POINTS, 0, m_centers.size() );
	}

}

// Draws a spiral
void GLViewer::draw()
{
	fastDraw();
}

void GLViewer::init()
{
  // Restore previous viewer state.
  restoreStateFromFile();
  update_buffers();
}

void GLViewer::update_buffers()
{
	const Particles &p = m_offline.particles() ;
	m_centers.reset( p.count(), p.centers().data(), GL_STATIC_DRAW )  ;

	// Colors
	Eigen::Matrix4Xf colors( 4, p.count() ) ;
	colors.topRows(2).setZero() ;
	colors.row(2) = p.masses().cast< float >() * 1./p.masses().maxCoeff() ;
	colors.row(3).setOnes() ;
	m_colors.reset( p.count(), colors.data(), GL_STATIC_DRAW )  ;

}

void GLViewer::set_frame(unsigned frame)
{
	if( m_offline.load_frame( frame ) ) {
		m_currentFrame = frame ;
	}
	update_buffers();
}

void GLViewer::keyPressEvent(QKeyEvent *e)
{
	switch (e->key())
	{
	case Qt::Key_I :
		next_frame();
		updateGL() ;
		break ;
	case Qt::Key_Q :
		QApplication::exit( 0 ) ;
	default:
		QGLViewer::keyPressEvent(e);
	}
}

} //d6
