

#ifndef D6_GLVIEWER_HH
#define D6_GLVIEWER_HH

#include "GLCamera.hh"

#include "visu/Offline.hh"

#include "gl/VertexBuffer.hh"
#include "gl/Shader.hh"
#include "gl/GrainRenderer.hh"
#include "gl/ShapeRenderer.hh"
#include "gl/Texture.hh"

#include <Eigen/Eigenvalues>

namespace d6
{
class Offline;

class GLViewer
{

public:
    GLViewer(const Offline &offline,
             const int width, const int height)
        : m_offline(offline), m_width(width), m_height(height)
    {
    }

    int width() const { return m_width; }
    int height() const { return m_height; }

    void init();

    void update_buffers();

    void draw() const;

    void snap();
    void frameAll();

    void rotate(float xAmount, float yAmount);
    void translate(float xAmount, float yAmount);
    void zoom(float amount);

private:
    const Offline &m_offline;
    
    int m_width;
    int m_height;

	bool 	 m_drawParticles  = true;
	bool 	 m_enableBending  = false;
	bool	 m_fastDraw  = false;
	bool 	 m_drawObjects = true;
	bool 	 m_drawOrientations  = false ;

    Camera m_camera;

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
};

} // namespace d6

#endif
