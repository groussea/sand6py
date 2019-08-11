#include "GLViewer.hh"

#include "visu/Offline.hh"

#include "geo/FieldBase.impl.hh"
#include "geo/Particles.hh"
#include "geo/Tensor.hh"
#include "geo/MeshImpl.hh"
#include "geo/LevelSet.hh"

#include "utils/Log.hh"

namespace d6
{

void GLViewer::init()
{
    GLint bufs, samples;
    glGetIntegerv(GL_SAMPLE_BUFFERS, &bufs);
    glGetIntegerv(GL_SAMPLES, &samples);
    Log::Debug() << "Using " << bufs << " buffers and " << samples << " samples" << std::endl;

    frameAll();
}

void GLViewer::update_buffers()
{
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
}

void GLViewer::draw() const
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_ACCUM_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

    {
        glEnableClientState(GL_VERTEX_ARRAY);
        m_centers.set_vertex_pointer(3);

        if (m_drawParticles)
        {
            gl::ColorPointer cp(m_colors);
            glDrawArrays(GL_POINTS, 0, m_centers.size());
        }
        else
        {
            glColor3f(0, 0, 1);
            glDrawArrays(GL_POINTS, 0, m_centers.size());
        }

        glDisableClientState(GL_VERTEX_ARRAY);
    }
}

void GLViewer::snap()
{
}


void GLViewer::frameAll()
{
    Log::Debug() << "Framing " << m_centers.size() << " particles" << std::endl;
    Eigen::Vector3f box = m_offline.box().cast<float>();

    Eigen::Vector3f position;
    position[0] = 0.5*box[0] ;
    position[1] = 3*box[1];
    position[2] = 0.5*box[2];

    m_camera.lookAt(position, 0.5*box, Eigen::Vector3f(0,0,1));
    m_camera.setPerspective(2*M_PI/3, m_width/(float)m_height, 0.01, 10*box.norm() );
    m_camera.apply();
}

void GLViewer::move(float xAmount, float yAmount)
{
    m_camera.rotate(-M_PI*xAmount/m_width, -M_PI*yAmount/m_height);
    m_camera.apply();
}

void GLViewer::zoom(float amount)
{
    m_camera.zoom(amount);
    m_camera.apply();
}

} // namespace d6