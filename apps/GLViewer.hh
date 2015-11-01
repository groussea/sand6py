#ifndef D6_GLVIEWER_HH
#define D6_GLVIEWER_HH

#include "visu/Offline.hh"
#include "gl/VertexBuffer.hh"
#include "gl/Shader.hh"

#include <QGLViewer/qglviewer.h>

namespace d6 {

class GLViewer : public QGLViewer
{

public:

	GLViewer( Offline & offline ) :
		m_offline( offline ), m_currentFrame(-1),
		m_enableBending( false )
	{

	}

	void set_frame( unsigned frame )  ;

	bool next_frame() {
		unsigned nextFrame = m_currentFrame + 1  ;
		set_frame( nextFrame ) ;
		return nextFrame == m_currentFrame ;
	}

protected :
  virtual void fastDraw();
  virtual void draw();
  virtual void init();
  virtual void animate();

  virtual void keyPressEvent(QKeyEvent *e);

private:
	void update_buffers() ;

	Offline& m_offline ;
	unsigned m_currentFrame ;

	bool m_enableBending ;

	gl::VertexBuffer3d m_centers ;
	gl::VertexBuffer4f m_colors  ;

	gl::VertexBuffer3f m_glyph ;
	gl::IndexBuffer	   m_glyphQuadIndices ;

	Eigen::Matrix< float, 16, Eigen::Dynamic> m_matrices ; // FIXME
	gl::VertexBuffer16f m_frames  ;

	Eigen::VectorXf m_densities ;
	gl::ArrayBufferf m_alpha ;

	Shader m_shader ;

} ;

} //d6

#endif
