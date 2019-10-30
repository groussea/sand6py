#include "GLViewer.hh"

#include "visu/Offline.hh"

#include "geo/FieldBase.impl.hh"
#include "geo/Particles.hh"
#include "geo/Tensor.hh"
#include "geo/MeshImpl.hh"
#include "geo/LevelSet.hh"

#include "utils/Log.hh"

#define fb_width  (2*2560)
#define fb_height (2*1440)

namespace d6
{

GLViewer::GLViewer(const Offline &offline,
                   const int nSamples,
                   const int width, const int height)
    : m_offline(offline), m_width(width), m_height(height),
      m_drawParticles( 0 == nSamples ), 
      m_grainsRenderer( offline, m_shapeRenderer, nSamples )
{
}

void GLViewer::init()
{
    GLint bufs, samples;
    glGetIntegerv(GL_SAMPLE_BUFFERS, &bufs);
    glGetIntegerv(GL_SAMPLES, &samples);
    Log::Debug() << "Using " << bufs << " buffers and " << samples << " samples" << std::endl;

    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

    frameAll();

	if( m_grainsRenderer.sampler().mode() == Sampler::VelocityCut ) {
        // TODO set orthographic
	} else {
	}

    config_shaders();

	m_shapeRenderer.init();

	if( renderSamples() ) {

		m_grainsRenderer.init();

		m_depthBuffer.reset( fb_width, fb_height );
		m_depthTexture.reset( GL_TEXTURE_2D );
		m_depthTexture.bind() ;

		float* data = 0 ;
		glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, m_depthBuffer.width(), m_depthBuffer.height(), 0,GL_DEPTH_COMPONENT, GL_FLOAT, data);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	}

	update_buffers();
}

void GLViewer::config_shaders()
{
	m_particlesShader.add_attribute("vertex") ;
	m_particlesShader.add_attribute("frame") ;
	m_particlesShader.add_attribute("density") ;
	m_particlesShader.add_uniform("model_view") ;
	m_particlesShader.add_uniform("projection") ;
	m_particlesShader.add_uniform("light_pos") ;
	m_particlesShader.load("particles_vertex", "particles_fragment") ;
	
	m_pointShader.add_attribute("vertex") ;
	m_pointShader.add_uniform("model_view") ;
	m_pointShader.add_uniform("projection") ;
    m_pointShader.load("point_vertex", "point_fragment") ;

	m_testShader.add_attribute("vertex") ;
	m_testShader.add_uniform("in_texture");
	//m_testShader.load("textest_vertex","textest_fragment") ;
}

void GLViewer::update_buffers()
{
    if(!m_camera.valid()) return;

    const Particles &p = m_offline.particles();
    m_centers.reset(p.count(), p.centers().data(), GL_STATIC_DRAW);

    if (m_drawParticles)
    {

        m_matrices.resize(16, p.count());
        m_densities.resize(p.count());

        // Compute movel-view matrix from tensor
        Eigen::Matrix4f mat;
        Mat tensor;

#pragma omp parallel for private(mat, tensor)
        for (size_t i = 0; i < p.count(); ++i)
        {
            mat.setIdentity();

            if (m_drawOrientations)
            {
                tensor_view(p.orient().col(i)).get(tensor);
                Eigen::SelfAdjointEigenSolver<Mat> es(tensor);

                const Vec ev = es.eigenvalues().array().max(0).sqrt();

                mat.block<3, 3>(0, 0) = (es.eigenvectors() * ev.asDiagonal()).cast<GLfloat>() * .5 * std::pow(p.volumes()[i], 1. / 3);
                mat.block<3, 1>(0, 3) = p.centers().col(i).cast<GLfloat>();

                m_densities[i] = p.volumes()[i];
            }
            else
            {
                tensor_view(p.frames().col(i)).get(tensor);
                Eigen::SelfAdjointEigenSolver<Mat> es(tensor);

                const Vec ev = es.eigenvalues().array().max(0).sqrt();
                const Scalar vol = 8 * ev.prod();

                mat.block<3, 3>(0, 0) = (es.eigenvectors() * ev.asDiagonal()).cast<GLfloat>();
                mat.block<3, 1>(0, 3) = p.centers().col(i).cast<GLfloat>();

                m_densities[i] = p.volumes()[i] / vol;
            }

            m_matrices.col(i) = Eigen::Matrix<GLfloat, 16, 1>::Map(mat.data(), mat.size());
        }
        m_frames.reset(p.count(), m_matrices.data(), GL_STATIC_DRAW);
        m_alpha.reset(p.count(), m_densities.data(), GL_STATIC_DRAW);

        // Colors
        Eigen::Matrix4Xf colors(4, p.count());
        colors.topRows(2).setZero();
        colors.row(2) = m_densities.cast<float>();
        colors.row(3).setOnes();
        m_colors.reset(p.count(), colors.data(), GL_STATIC_DRAW);
    }

	if( renderSamples() ) {
		m_grainsRenderer.update_buffers();
	}

    update_vaos();
}

void GLViewer::update_vaos()
{

    if(m_pointShader.ok()){
        gl::ArrayObject::Using vao(m_pointArrays);
        //Vertices
        gl::VertexAttribPointer vap(m_centers, m_pointShader.attribute("vertex"));
    }

    if(m_drawParticles) {
        if(m_particlesShader.ok()){
            gl::ArrayObject::Using vao(m_shapeRenderer.sphereVertexArrays());
           
            //Vertices
            gl::VertexAttribPointer vap(m_shapeRenderer.sphereVertices(), m_particlesShader.attribute("vertex"));
            // Densities
            gl::VertexAttribPointer ap(m_alpha, m_particlesShader.attribute("density"), false, 1);
            //Frames
            gl::ArrayAttribPointer<4> fp(m_frames, m_particlesShader.attribute("frame"), false, 1);
        }
    }
	
    m_shapeRenderer.clear_buffers();
    for( const LevelSet::Ptr& ls: m_offline.levelSets() ) {
        m_shapeRenderer.setup_buffers(*ls, m_offline.box().cast<float>());
    }
}

Eigen::Vector3f GLViewer::lightPosition() const
{
	return ( m_offline.box() / 2 ).cast< float >() + 3 * m_offline.box().norm() * m_lightDirection.normalized() ;
}

void GLViewer::draw() const
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT |  GL_STENCIL_BUFFER_BIT);

    if(m_fastDraw)
    {
        if( m_pointShader.ok() ) {
            glPointSize(2.0f);
			UsingShader sh( m_pointShader ) ;

			// Model-view
            sh.bindMVP(m_camera.viewMatrix.data(), m_camera.projectionMatrix.data());

			//Vertices
            gl::ArrayObject::Using vao(m_pointArrays);
            glDrawArrays( GL_POINTS, 0, m_centers.size());

		}  else {
            std::cerr << "Invalid point shader" << std::endl; 
        }
        return;
    }

    glEnable(GL_DEPTH_TEST);

	if( m_enableBending ) {
		glEnable (GL_BLEND);
		glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}

	bool shadowed = false ;

	if( m_drawParticles )
	{
		 if( m_particlesShader.ok() ) {

			UsingShader sh( m_particlesShader ) ;

            Eigen::Vector3f lightPos = lightPosition();
            glUniform3fv(m_particlesShader.uniform("light_pos"), 1, lightPos.data());

			// Model-view
            sh.bindMVP(m_camera.viewMatrix.data(), m_camera.projectionMatrix.data());
            
            gl::ArrayObject::Using vao(m_shapeRenderer.sphereVertexArrays());
			glDrawElementsInstanced( GL_TRIANGLES, m_shapeRenderer.sphereTriIndices().size(), GL_UNSIGNED_INT, 0, m_matrices.cols() );

		}  else {
            std::cerr << "Invalid particles shader" << std::endl; 
        }

	}

    Camera depthCam = m_camera;
    Camera viewCam = m_camera;

    if( renderSamples() )
	{

		// Set-up light POV camera
        const Eigen::Vector3f &light_pos = lightPosition();
        Eigen::Vector3f box = m_offline.box().cast<float>();

        const Eigen::Vector3f sceneCenter = 0.5f*box;
        const float centerDepth = (light_pos - sceneCenter).norm(); 
		depthCam.lookAt(light_pos, sceneCenter, depthCam.up);
        float fov = std::atan2(box.norm(), centerDepth);
        depthCam.setPerspective(fov, fb_width/(float)fb_height, 0.5*centerDepth, 2*centerDepth);

        const Eigen::Matrix4f depthMVP = depthCam.projectionMatrix * depthCam.viewMatrix;


        if(  m_grainsRenderer.sampler().mode() != Sampler::VelocityCut )
		{

			// Compute grains shadowing
			UsingFrameBuffer fb( m_depthBuffer ) ;

			glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, m_depthTexture.id(), 0);
			glDrawBuffer(GL_NONE); // No color buffer is drawn to.

			if( m_depthBuffer.check_complete() )
				Log::Error() << "Frame buffer incomplete" << std::endl ;

			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

			const Scalar pixelSize = false /*cam.type() == qglviewer::Camera::ORTHOGRAPHIC*/
			        ? 0
			        : fb_height / std::tan(depthCam.fieldOfView / 2)  ;

			m_grainsRenderer.compute_shadow( pixelSize, depthMVP );

			if( m_drawObjects ) {
				for( const LevelSet::Ptr& ls: m_offline.levelSets() ) {
					m_shapeRenderer.compute_shadow( *ls, depthCam.viewMatrix, depthCam.projectionMatrix ) ;
				}
			}

			shadowed = true ;
		}

		if(0){
			UsingShader sh( m_testShader ) ;

			UsingTexture tx( m_depthTexture ) ;
			tx.bindUniform( m_testShader.uniform("in_texture") );

			gl::VertexAttribPointer vap_v( m_shapeRenderer.squareVertices(), m_testShader.attribute("vertex") ) ;
			gl::VertexPointer vp( m_shapeRenderer.squareVertices() ) ;
			glDrawArrays( GL_TRIANGLES, 0, m_shapeRenderer.squareVertices().size() ) ;

		}

		if(1){
			const Scalar pixelSize = false /*cam.type() == qglviewer::Camera::ORTHOGRAPHIC*/
			        ? 0
			        : m_height / std::tan(viewCam.fieldOfView / 2)  ;

			m_grainsRenderer.draw( m_depthTexture, lightPosition(), viewCam.viewMatrix, viewCam.projectionMatrix, pixelSize, depthMVP  );
		}

		glDisable( GL_PROGRAM_POINT_SIZE ) ;
	}

	glDisable (GL_BLEND);

	if(m_drawObjects) {
		for( const LevelSet::Ptr& ls: m_offline.levelSets() ) {
			m_shapeRenderer.draw( *ls, m_offline.box(), lightPosition(), shadowed, m_depthTexture, 
            viewCam.viewMatrix, viewCam.projectionMatrix, depthCam.viewMatrix, depthCam.projectionMatrix );
		}
	}

}

void GLViewer::snap()
{
}


void GLViewer::frameAll()
{
    Log::Debug() << "Framing " << m_centers.size() << " particles" << std::endl;
    Eigen::Vector3f box = m_offline.box().cast<float>();

    if (!m_camera.valid())
    {
        Eigen::Vector3f position;
        position[0] = 0.5 * box[0];
        position[1] = 3 * box[1];
        position[2] = 0.5 * box[2];
        m_camera.lookAt(position, 0.5 * box, Eigen::Vector3f(0, 0, 1));
        m_camera.setPerspective(M_PI/3, m_width / (float)m_height, 0.1*box.norm(), 10 * box.norm());
    }
    else
    {
        m_camera.frame(box);
    }

    m_camera.apply();
}

void GLViewer::rotate(float xAmount, float yAmount)
{
    m_camera.rotate(-M_PI*xAmount/m_width, -M_PI*yAmount/m_height);
    m_camera.apply();
}

void GLViewer::translate(float xAmount, float yAmount)
{
    m_camera.translate(xAmount/m_width, -yAmount/m_height);
    m_camera.apply();
}

void GLViewer::zoom(float amount)
{
    m_camera.zoom(amount);
    m_camera.apply();
}

} // namespace d6