#ifndef D6_GLVIEWER_HH
#define D6_GLVIEWER_HH

#include "visu/Offline.hh"
#include "visu/VertexBuffer.hh"

#include <QGLViewer/qglviewer.h>

namespace d6 {

class GLViewer : public QGLViewer
{

public:

	GLViewer( Offline & offline ) :
		m_offline( offline ), m_currentFrame(-1)
	{

	}

	void set_frame( unsigned frame )  ;

	void next_frame() {
		set_frame( m_currentFrame + 1 ) ;
	}

protected :
  virtual void fastDraw();
  virtual void draw();
  virtual void init();

  virtual void keyPressEvent(QKeyEvent *e);

private:
	void update_buffers() ;

	Offline& m_offline ;
	unsigned m_currentFrame ;

	gl::VertexBuffer3d m_centers ;
	gl::VertexBuffer4f m_colors  ;

} ;

} //d6

#endif
