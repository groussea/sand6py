#ifndef D6_GLVIEWER_HH
#define D6_GLVIEWER_HH

#include "visu/Offline.hh"
#include "visu/Sampler.hh"

#include "gl/VertexBuffer.hh"
#include "gl/Shader.hh"

#include <QGLViewer/qglviewer.h>

namespace d6 {

class GLViewer : public QGLViewer
{

public:

	GLViewer( Offline & offline, bool sample = true ) :
		m_offline( offline ), m_sampler( offline ), m_renderSamples( sample ),
		m_currentFrame(-1),
		m_drawParticles( !sample ), m_enableBending( false ), m_snapshotting(false),
		m_lastSnapped( m_currentFrame )
	{

	}

	void set_frame( unsigned frame )  ;

	bool next_frame() {
		unsigned nextFrame = m_currentFrame + 1  ;
		set_frame( nextFrame ) ;
		return nextFrame == m_currentFrame ;
	}
	bool prev_frame() {
		unsigned nextFrame = m_currentFrame - 1  ;
		if( m_currentFrame > 0 )
			set_frame( nextFrame ) ;
		return nextFrame == m_currentFrame ;
	}

protected :
  virtual void fastDraw();
  virtual void draw();
  virtual void drawWithNames();
  virtual void init();
  virtual void animate();

  virtual void keyPressEvent(QKeyEvent *e);
  virtual void postSelection(const QPoint& ) ;

private:

	void update_buffers() ;
	void drawObject( const LevelSet& ls ) ;
	void snap() ;

	Offline& m_offline ;
	Sampler  m_sampler ;

	bool  	 m_renderSamples  ;
	unsigned m_currentFrame ;

	bool 	 m_drawParticles ;
	bool 	 m_enableBending ;
	bool 	 m_snapshotting ;
	unsigned m_lastSnapped ;

	gl::VertexBuffer3d m_centers ;
	gl::VertexBuffer4f m_colors  ;

	gl::VertexBuffer3f m_glyph ;
	gl::IndexBuffer	   m_glyphQuadIndices ;

	Eigen::Matrix< float, 16, Eigen::Dynamic> m_matrices ; // FIXME
	gl::VertexBuffer16f m_frames  ;

	Eigen::VectorXf m_densities ;
	gl::ArrayBufferf m_alpha ;
	
	gl::VertexBuffer3f m_grainVertices ;
	gl::VertexBuffer3f m_grainNormals ;
	gl::ArrayBufferf m_grainVisibility ;

	Shader m_shader ;
	Shader m_grainsShader ;

} ;

} //d6

#endif
