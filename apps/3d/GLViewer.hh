

#ifndef D6_GLVIEWER_HH
#define D6_GLVIEWER_HH

#include "GLCamera.hh"

#include "visu/Offline.hh"

#include "gl/VertexBuffer.hh"
#include "gl/Shader.hh"
#include "gl/GrainRenderer.hh"
#include "gl/ShapeRenderer.hh"
#include "gl/Texture.hh"
#include "gl/gImage.hh"

#include <Eigen/Eigenvalues>
#include "geo/LevelSet.hh"
#include "geo/LevelSet.io.hh"

#include "visu/Axes.hh"

#include <memory>

namespace d6
{
class Offline;
class LevelSet ;

class GLViewer
{

public:
    
    explicit GLViewer(const Offline &offline,
             const int nSamples,
             const int width, const int height);
        GLViewer(const GLViewer&) = delete;
        GLViewer& operator=(const GLViewer&) = delete;
        ~GLViewer() = default;


    int width() const { return m_width; }
    int height() const { return m_height; }

    void init();
    void config_shaders();

    void update_buffers();
    void update_vaos();

    void draw() const;

    void snap(unsigned  m_currentFrame);
    void frameAll();

    void rotate(float xAmount, float yAmount);
    void translate(float xAmount, float yAmount);
    void zoom(float amount);
    void look_at(Eigen::Vector3f lookPos);
    void set_cam_pos(Eigen::Vector3f camPos);
    Eigen::Vector3f  get_cam_pos() const;
    void setFastDraw(bool enable)
    {
        m_fastDraw = enable;
    }
    void toggleFastDraw() {
        m_fastDraw = !m_fastDraw;
    }

	bool renderSamples() const {
		return m_grainsRenderer.valid() ;
	}

	GrainRenderer& grainsRenderer() { return m_grainsRenderer ;}




private:
    Eigen::Vector3f lightPosition() const;

    const Offline &m_offline;

    const LevelSet& levelSet() const
	{ return *m_levelSet ; }
	const LevelSet* levelSetPtr() const
	{ return m_levelSet.get() ; }

    std::unique_ptr< LevelSet > m_levelSet ;

	std::vector< Axes >  m_axes ;
    // std::unique_ptr< LevelSet >  m_axes;

    int m_width;
    int m_height;

	bool 	 m_drawParticles  = true;
	bool 	 m_enableBending  = true;
	bool	 m_fastDraw  = false;
	bool 	 m_drawObjects = true;
	bool 	 m_drawOrientations  = false ;
    bool     m_drawAxis = true;

    Eigen::Vector3f m_lightDirection = Eigen::Vector3f(0.5, 0.5, 1);

    Camera m_camera;

	gl::VertexBuffer3d m_centers ;
	gl::VertexBuffer4f m_colors  ;
    gl::ArrayObject m_pointArrays;

	Eigen::Matrix< float, 16, Eigen::Dynamic> m_matrices ;
	gl::VertexBuffer16f m_frames  ;

	Eigen::VectorXf m_densities ;
	gl::ArrayBufferf m_alpha ;

	Shader m_particlesShader ;
	Shader m_pointShader ;
	Shader m_testShader ;

	Texture     m_depthTexture ;
	FrameBuffer m_depthBuffer  ;

	ShapeRenderer m_shapeRenderer ;
	GrainRenderer m_grainsRenderer ;
};

} // namespace d6

#endif
