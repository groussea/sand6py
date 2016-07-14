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

#ifndef D6_GLVIEWER_HH
#define D6_GLVIEWER_HH

#include "visu/Offline.hh"

#include "gl/VertexBuffer.hh"
#include "gl/Shader.hh"
#include "gl/GrainRenderer.hh"
#include "gl/ShapeRenderer.hh"
#include "gl/Texture.hh"

#include <QGLViewer/qglviewer.h>

namespace d6 {

class GLViewer : public QGLViewer
{

public:

	GLViewer( const QGLFormat& glFormat, Offline & offline, unsigned nSamples,
			unsigned width, unsigned height
		) ;

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

	bool renderSamples() const {
		return m_grainsRenderer.valid() ;
	}

	GrainRenderer& grainsRenderer() {
		return m_grainsRenderer ;
	}
	const GrainRenderer& grainsRenderer() const {
		return m_grainsRenderer ;
	}

		void setLightDirection( const Eigen::Vector3f& dir )
		{
		  m_lightDirection = dir ;
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

	Eigen::Vector3f lightPosition() const ;

	void update_buffers() ;
	void drawObject( const LevelSet& ls ) ;
	void snap() ;

	Offline& m_offline ;

		unsigned m_vp_width ;
		unsigned m_vp_height ;
		Eigen::Vector3f m_lightDirection ;

	unsigned m_currentFrame ;

	bool 	 m_drawParticles ;
	bool 	 m_enableBending ;
	bool	 m_fastDraw ;
	bool 	 m_drawObjects ;
	bool 	 m_drawOrientations ;
	bool 	 m_snapshotting ;

	unsigned m_lastSnapped ;

	gl::VertexBuffer3d m_centers ;
	gl::VertexBuffer4f m_colors  ;

	Eigen::Matrix< float, 16, Eigen::Dynamic> m_matrices ;
	gl::VertexBuffer16f m_frames  ;

	Eigen::VectorXf m_densities ;
	gl::ArrayBufferf m_alpha ;

	Shader m_particlesShader ;
	Shader m_testShader ;

	Texture     m_depthTexture ;
	FrameBuffer m_depthBuffer  ;

	ShapeRenderer m_shapeRenderer ;
	GrainRenderer m_grainsRenderer ;

} ;

} //d6

#endif
