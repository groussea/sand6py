/*
 * This file is part of Sand6, a C++ continuum-based granular simulator.
 *
 * Copyright 2016 Gilles Daviet <gilles.daviet@inria.fr> (Inria - Universit√© Grenoble Alpes)
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

#include "GLViewer.hh"

#include "geo/FieldBase.impl.hh"
#include "geo/MeshImpl.hh"
#include "geo/LevelSet.hh"

#include "mono/Phase.hh"
#include "diphasic/FluidPhase.hh"
#include "visu/Offline.hh"

#include "utils/Log.hh"

#include <Eigen/Eigenvalues>

#include <iostream>
#include <cstdio>

#define TEX_S 512u

namespace d6 {

void GLViewer::move(double x, double y)
{
	m_xOffset -= x / m_width ;
	m_yOffset += y / m_height ;
	updateViewport();
}

void GLViewer::zoom(double factor)
{
	m_zoom *= factor ;
	m_xOffset = ( m_xOffset + .5 )/factor - .5 ;
	m_yOffset = ( m_yOffset + .5 )/factor - .5 ;
	updateViewport();
}
#include "geo/FieldBase.impl.hh"


void GLViewer::updateViewport()
{
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();

	const Config& c = m_offline.config() ;

	const double relWidth  = m_zoom * c.box[0] ; //*m_width)/m_height ;
	const double relHeight = m_zoom * c.box[1] ;

	glOrtho ( (m_xOffset-.1)*relWidth, (m_xOffset+1.1)*relWidth, (m_yOffset-.1)*relHeight, (m_yOffset+1.1)*relHeight, 0, 1);
	glMatrixMode (GL_MODELVIEW);
}


template < typename Scalar >
static Scalar clamp( const Scalar v, const Scalar min=0, const Scalar max=1) {
	return std::min( max, std::max( min, v )) ;
}

inline float rnd( float* data, float x, float y )
{

	const unsigned px = clamp( (unsigned) std::floor( x ), 0u, TEX_S-1) ;
	const unsigned py = clamp( (unsigned) std::floor( y ), 0u, TEX_S-1) ;

	const unsigned off = px + py*TEX_S ;

	x -= px ;
	y -= py ;

	return  data[ off ] * (1-x) * (1-y)  +
	        data[ 1 + off ] * y * (1-x) +
	        data[ off + TEX_S ] * (1-y) * x +
	        data[ 1 + off + TEX_S ] * x * y ;
}

void GLViewer::snap()
{
	using namespace std ;

	const unsigned g_gl_width = m_width ;
	const unsigned g_gl_height = m_height ;
	//http://antongerdelan.net/opengl/screencapture.html
	unsigned char* buffer = (unsigned char*)malloc (g_gl_width * g_gl_height * 3);
	glReadPixels (0, 0, g_gl_width, g_gl_height, GL_RGB, GL_UNSIGNED_BYTE, buffer) ;
	char name[1024];
	unsigned t = m_snapId++ ;

	if( m_snapId > 9999 ) std::exit(1) ;

	sprintf (name, "snapshot-%04d.raw", t);
	std::cout << "Snapping " << name << std::endl ;
	FILE* fp = fopen (name, "wb");
	if (!fp) {
		fprintf (stderr, "ERROR: could not open %s for writing\n", name);
		return ;
	}
	int bytes_in_row = g_gl_width * 3;
	int bytes_left = g_gl_width * g_gl_height * 3;
	while (bytes_left > 0) {
		int start_of_row = bytes_left - bytes_in_row;
		fwrite (&buffer[start_of_row], 1, bytes_in_row, fp);
		bytes_left -= bytes_in_row;
	}
	fclose (fp);
	free (buffer);
	char command[2048];
	sprintf (
	            command,
	            "convert -depth 8 -size %ix%i rgb:snapshot-%04d.raw snapshot-%04d.png",
	            g_gl_width,
	            g_gl_height,
	            t,
	            t
	            );
	system (command);
	sprintf (command, "rm -f %s", name);
	system (command);
}


void GLViewer::init( )
{
//	glEnable(GL_MULTISAMPLE);
	GLint bufs, samples ;
	glGetIntegerv( GL_SAMPLE_BUFFERS, &bufs ) ;
	glGetIntegerv( GL_SAMPLES, &samples ) ;
	Log::Debug() << "Using " << bufs << " buffers and " << samples << " samples" << std::endl ;

	updateViewport();

	glClearColor(1.,1.,1.,1) ;

	// Texture
	{
		m_rndData.resize( TEX_S * TEX_S  );
		m_texData.resize( TEX_S * TEX_S * 3 );
		Eigen::Matrix< float, Eigen::Dynamic, 1 >
		        ::Map( &m_rndData[0], TEX_S * TEX_S )
		        = ( Eigen::VectorXf::Random( TEX_S * TEX_S ).array().abs().matrix() )  ;

		glGenTextures(1, &m_texId);
		glBindTexture( GL_TEXTURE_2D, m_texId );
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	}

	//

	m_shouldRender[ eGrid ] = true ;
	m_shouldRender[ eColors ] = true ;
	m_shouldRender[ eVectors ] = true ;
	m_shouldRender[ eTensors ] = false ;
	m_shouldRender[ eStream ] = false ;
	m_shouldRender[ eParticles ] = true ;

	m_particleShader.add_uniform("model_view");
	m_particleShader.add_uniform("projection");
	m_particleShader.add_attribute("vertex");
	m_particleShader.add_attribute("density");
	m_particleShader.add_attribute("value");
	m_particleShader.add_attribute("frame");
	m_particleShader.load("2d/particle_vertex","2d/particle_fragment") ;

	m_vectorShader.add_uniform("model_view");
	m_vectorShader.add_uniform("projection");
	m_vectorShader.add_attribute("vertex");
	m_vectorShader.add_attribute("value");
	m_vectorShader.load("2d/vector_vertex","2d/vector_fragment") ;

	m_tensorShader.add_uniform("model_view");
	m_tensorShader.add_uniform("projection");
	m_tensorShader.add_attribute("vertex");
	m_tensorShader.add_attribute("value");
	m_tensorShader.load("2d/tensor_vertex","2d/tensor_fragment") ;
}

void GLViewer::update_buffers()
{
	const MeshType& g = m_offline.meshes().primal() ;
	const typename VectorField::ShapeFuncImpl shape( g ) ;

	// Grid nodes
	// Grid quad indices
	{
		DynMatW vtx ;
		DynMati ids ;
		shape.derived().build_visu_mesh( vtx, ids ) ;

		Eigen::Matrix2Xf vertices = vtx.cast< float >() ;

		const Index NV     = ids.rows() ;
		const Index nCells = ids.cols() ;

		std::vector<unsigned> nodeIndices( NV * nCells ) ;

		for( Index i = 0 ; i < nCells ; ++i )
		{
			if( NV == 4 ) {
				nodeIndices[ 4*i + 0 ] = ids(0,i) ;
				nodeIndices[ 4*i + 1 ] = ids(1,i) ;
				nodeIndices[ 4*i + 2 ] = ids(3,i) ;
				nodeIndices[ 4*i + 3 ] = ids(2,i) ;
			} else {
				for( int k = 0 ; k < NV ; ++k ) {
					nodeIndices[ NV*i + k ] = ids(k,i) ;
				}
			}
		}
		m_gridPrimitive = NV == 4 ? GL_QUADS : GL_TRIANGLES ;

		m_gridVertices.reset( vertices.cols(), vertices.data(), GL_STATIC_DRAW  ) ;
		m_gridQuadIndices.reset( nodeIndices.size(), &nodeIndices[0], GL_STATIC_DRAW  ) ;
	}


	if( m_shouldRender[ eParticles ])
		update_particle_buffers();

	if( m_shouldRender[ eColors ])
		update_color_buffers();

	if( m_shouldRender[ eVectors ])
		update_vector_buffers();

	if( m_shouldRender[ eTensors ])
		update_tensor_buffers();

	if( m_shouldRender[ eStream ])
		update_texture();

}

void GLViewer::update_particle_buffers()
{
	//Particles
	const Particles &p = m_offline.particles() ;
	m_particles.reset( p.count(), p.centers().data(), GL_STATIC_DRAW )  ;

	//		m_matrices.resize( 9, p.count() );
	Eigen::VectorXf densities ( p.count() );
	Eigen::VectorXf velocities ( p.count() );
	Eigen::Matrix< float, 4, Eigen::Dynamic > frames (4, p.count()) ;

	// Compute movel-view matrix from tensor
	Mat tensor ;
	Eigen::Matrix2f mat ;

#pragma omp parallel for private(tensor, mat)
	for( size_t i = 0 ; i < p.count() ; ++i ) {

		if( m_particleShape == psOrientation) {
			tensor_view( p.orient().col( i ) ).get( tensor ) ;
			Eigen::SelfAdjointEigenSolver<Mat> es( tensor );

			const Vec ev = es.eigenvalues().array().max(0).sqrt().max( 1.e-1 ) ;

			mat = ( es.eigenvectors() * ev.asDiagonal() ).cast< GLfloat >()
			        * std::pow( p.volumes()[i], 1./2 )  ;

			densities[i] = p.volumes()[i]  ;

		} else {
			tensor_view( p.frames().col( i ) ).get( tensor ) ;
			Eigen::SelfAdjointEigenSolver<Mat> es( tensor );

			const Vec ev = es.eigenvalues().array().max(0).sqrt().max( 1.e-2 ) ;
			const Scalar vol = 4 * ev.prod() ;

			if( m_particleShape == psFrame )
				mat = ( es.eigenvectors() * ev.asDiagonal() ).cast< GLfloat >()  ;
			else
				mat = std::sqrt(p.volumes()[i]/M_PI) * Eigen::Matrix2f::Identity() ;

			densities[i] = p.volumes()[i] / vol ;
		}

		velocities[i] = p.velocities().col(i).norm() ;


		frames.col(i) = Eigen::Matrix< GLfloat, 4, 1 >::Map( mat.data() ) ;

	}

	m_particleFrames.reset( p.count(), frames.data(), GL_STATIC_DRAW )  ;
	m_particleAlpha.reset ( p.count(), densities.data(), GL_STATIC_DRAW )  ;
	m_particleColors.reset ( p.count(), velocities.data(), GL_STATIC_DRAW )  ;

	Log::Verbose() << "Particle max vel \t" << velocities.maxCoeff() << std::endl ;
}

// Grid colors
void GLViewer::update_color_buffers()
{
	const auto & field = getScalarEntity() ;

	Eigen::Matrix<float, 3, Eigen::Dynamic > colors( 3, field.size() ) ;
	colors.setZero() ;

	if( m_scalarEntity == seFluidVel ) {

#pragma omp parallel for
		for( Index i = 0 ; i < field.size() ; ++i ) {
			const Scalar vis = std::log( 1. + 100 * field[i] ) ;

			float mat_0 = clamp( 1 - 2*std::fabs(vis/5-1.0), 0.0, 1.0 ) ;
			float mat_1 = clamp( 1 - 2*std::fabs(vis/5-0.5), 0.0, 1.0 ) ;
			float mat_2 = clamp( 1 - 2*std::fabs(vis/5-0.0), 0.0, 1.0 ) ;
			float mat_3 = clamp( 1 + 2*   (vis/5-1.5), 0.0, 1.0 ) ;

			colors.col(i) += .5*mat_0*Eigen::Vector3f(1.0, 0.0, 0.0);
			colors.col(i) += .5*mat_1*Eigen::Vector3f(0.0, 1.0, 0.0);
			colors.col(i) += .5*mat_2*Eigen::Vector3f(0.0, 0.0, 1.0);
			colors.col(i) += .5*mat_3*Eigen::Vector3f(1.0, 1.0, 1.0);
		}

	} else {
		const double max_val = field.flatten().maxCoeff() ;
		const double min_val = field.flatten().minCoeff() ;
		if( max_val - min_val > 1.e-12 )
		{
			colors.row(1) = ( (field.flatten().array() - min_val)/(max_val - min_val) ).cast<float>() ;
		}

		Log::Verbose() << "Scalar range\t" << min_val << " ;\t" << max_val << std::endl ;
	}

	m_gridColors.reset( colors.cols(), colors.data() ) ;

}


// Arrows
void GLViewer::update_vector_buffers()
{
	const auto & field = getVectorEntity() ;

	Eigen::Matrix<float, 2, Eigen::Dynamic > vectors( 2, field.size() ) ;

	vectors.setZero() ;

	ScalarField norm = field.norm() ;

	const double max_val = norm.max_abs() ;
	if( max_val > 1.e-12 )
	{
		Eigen::VectorXf::Map( vectors.data(), vectors.size() ) = ( field.flatten().array() / max_val ).cast<float>() ;
	}

	m_arrows.reset( vectors.cols(), vectors.data() ) ;

	Log::Verbose() << "Vector max norm\t" <<  max_val << std::endl ;
}

// Tensors
void GLViewer::update_tensor_buffers()
{
	const auto & field = getTensorEntity() ;

	Eigen::Matrix<float, 3, Eigen::Dynamic > tensors( 3, field.size()) ;

	tensors.setZero() ;

	ScalarField norm = field.norm() ;

	const double max_val = norm.max_abs() ;
	if( max_val > 1.e-12 )
	{
		Eigen::VectorXf::Map( tensors.data(), tensors.size() ) = ( field.flatten().array() / max_val ).cast<float>() ;
	}

	m_tensors.reset( tensors.cols(), tensors.data() ) ;

	Log::Verbose() << "Tensor max norm\t" <<  max_val << std::endl ;
}

void GLViewer::update_texture()
{
	const Config& c = m_offline.config() ;
	const MeshType& g = m_offline.meshes().primal() ;

	// Texture

	glBindTexture( GL_TEXTURE_2D, m_texId );
	const float sx = c.box[0] /TEX_S ;
	const float sy = c.box[1] /TEX_S ;

	const VectorField& field = getVectorEntity() ;
	const double max_val = field.max_abs() ;

#pragma omp parallel for
	for( unsigned c = 0 ; c < TEX_S ; ++c )
	{
		for( unsigned r = 0 ; r < TEX_S ; ++r )
		{


			float x = c ;
			float y = r ;
			float v = rnd( &m_rndData[0], x, y ) ;

			unsigned kmax = 5 ;

			Vec u = field(g.locate(Vec(x*sx,y*sy))) ;
			float un = u.norm() ;

			float ax =  u[0] / max_val ;
			float ay =  u[1] / max_val ;
			float ampl = un / max_val ;

			for( unsigned k = 0 ; k < kmax ; ++k )
			{
				un = u.norm() ;
				if( un < 1.e-6 ) break ;

				x -= u[0] / un ;
				y -= u[1] / un ;

				v +=  rnd( &m_rndData[0], x, y ) ;
				u = field(g.locate(Vec(x*sx,y*sy))) ;
			}
			x = c ; y = r ;
			for( unsigned k = 0 ; k < kmax ; ++k )
			{
				u = field(g.locate(Vec(x*sx,y*sy))) ;
				un = u.norm() ;
				if( un < 1.e-6 ) break ;

				x += u[0] / un ;
				y += u[1] / un ;

				v +=  rnd( &m_rndData[0], x, y ) ;
			}
			v /= 2 * kmax ;

			const std::size_t idx = 3*( c + r*TEX_S );

			m_texData[idx]   = v*(1-ay) ;
			m_texData[idx+1] = v*(1-ax) ;
			m_texData[idx+2] = v*(1-ampl) ;

		}
	}
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, TEX_S, TEX_S, 0, GL_RGB, GL_FLOAT, &m_texData[0]);
}

void GLViewer::draw( ) const
{
	if(!m_gridQuadIndices.valid() ) return ;

	const Config& c = m_offline.config() ;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_ACCUM_BUFFER_BIT | GL_STENCIL_BUFFER_BIT );
	glLineWidth( 1.f );

	//BG
	glColor3f(0.8,0.8,0.8) ;
	glBegin(GL_QUADS);
	glVertex2f(0, 0);
	glVertex2f(0, c.box[1]);
	glVertex2f(c.box[0], c.box[1]);
	glVertex2f(c.box[0], 0);
	glEnd();

	// Grid

	if( m_shouldRender[ eColors ])
	{
		m_gridQuadIndices.bind() ;
		gl::VertexPointer vp( m_gridVertices ) ;
		gl::ColorPointer cp( m_gridColors ) ;
		glDrawElements( m_gridPrimitive, m_gridQuadIndices.size(), GL_UNSIGNED_INT, 0 );
	}

	glColor3f( .5, .5, .5 );
	if( m_shouldRender[ eGrid ])
	{
		m_gridQuadIndices.bind() ;
		gl::VertexPointer vp( m_gridVertices ) ;
		glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
		glDrawElements( m_gridPrimitive, m_gridQuadIndices.size(), GL_UNSIGNED_INT, 0 );
		glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

//		glDrawArrays( GL_LINES, 0, m_gridVertices.size()/2 );
		// Axis
		glLineWidth( 2.5f );
		glColor3f( 0., 1., 0. );
		glBegin( GL_LINES ) ; glVertex2f( 0, 0 ) ; glVertex2f( .1*c.units().fromSI(Units::Length), 0 ) ; glEnd() ;
		glColor3f( 0., 0., 1. );
		glBegin( GL_LINES ) ; glVertex2f( 0, 0 ) ; glVertex2f( 0, .1*c.units().fromSI(Units::Length) ) ; glEnd() ;
		glLineWidth( 1.f );
	}

	//Particles
	if(m_shouldRender[eParticles])
	{

		if( m_particleShader.ok() ) {

			UsingShader sh( m_particleShader ) ;
			sh.bindMVP();

			gl::VertexAttribPointer   vap( m_particles, m_particleShader.attribute("vertex"), false, 1) ;
			gl::VertexAttribPointer   dap( m_particleAlpha, m_particleShader.attribute("density"), false, 1) ;
			gl::VertexAttribPointer   cap( m_particleColors, m_particleShader.attribute("value"), false, 1) ;
			gl::VertexAttribPointer   aap( m_particleFrames, m_particleShader.attribute("frame"), false, 1) ;

			glDrawArraysInstanced( GL_QUADS, 0, 4, m_particles.size() );

		} else {
			glColor3f( 1., 0., 0. );
			glPointSize( 2.f ) ;

			gl::VertexPointer vp( m_particles ) ;
			glDrawArrays( GL_POINTS, 0, m_particles.size() );
		}
	}

	//Arrows
	if(m_shouldRender[eVectors])
	{
		if( m_vectorShader.ok() ) {

			UsingShader sh( m_vectorShader ) ;
			sh.bindMVP();

			gl::VertexAttribPointer   vap( m_gridVertices, m_vectorShader.attribute("vertex"), false, 1) ;
			gl::VertexAttribPointer   aap( m_arrows, m_vectorShader.attribute("value"), false, 1) ;

			glDrawArraysInstanced( GL_TRIANGLES, 0, 3, m_gridVertices.size() );
		}
		glEnableVertexAttribArray( 0 ); // ???
	}

	//Tensors
	if(m_shouldRender[eTensors])
	{
		if( m_tensorShader.ok() ) {

			UsingShader sh( m_tensorShader ) ;
			sh.bindMVP();

			gl::VertexAttribPointer  vap( m_gridVertices, m_tensorShader.attribute("vertex"), false, 1) ;
			gl::VertexAttribPointer  aap( m_tensors, m_tensorShader.attribute("value"), false, 1) ;

			glDrawArraysInstanced( GL_QUADS, 0, 4, m_gridVertices.size() );
		}
		glEnableVertexAttribArray( 0 ); // ???
	}

	// Texture
	if( m_shouldRender[ eStream ])
	{

		glBindTexture( GL_TEXTURE_2D, m_texId );

		glEnable(GL_TEXTURE_2D);
		glColor3f(1,1,1);
		glBegin(GL_QUADS);
		glTexCoord2f(0, 0); glVertex2f(0, 0);
		glTexCoord2f(0, 1); glVertex2f(0, c.box[1]);
		glTexCoord2f(1, 1); glVertex2f(c.box[0], c.box[1]);
		glTexCoord2f(1, 0); glVertex2f(c.box[0], 0);
		glEnd();
		glDisable(GL_TEXTURE_2D);
	}


	//LS
	for( const LevelSet::Ptr& ls: m_offline.levelSets() ) {
		m_shapeRenderer.draw( *ls, m_offline.box());
	}
}

typename GLViewer::VectorField GLViewer::getVectorEntity() const
{
	switch ( m_vectorEntity ) {
	case veVelocity :
		return m_offline.grains().velocity ;
	case veForce :
		return m_offline.grains().fcontact ;
	case veProj :
		return m_offline.grains().geo_proj ;
	case veFluidVel :
		return m_offline.fluid() .velocity ;
	case veFluctu :
		return GLViewer::VectorField(m_offline.grains().velocity).axpy( m_offline.fluid().velocity, -1 ) ;
	default:
		return m_offline.grains().grad_phi ;
	}
}
typename GLViewer::ScalarField GLViewer::getScalarEntity() const
{
	switch ( m_scalarEntity ) {
	case seFraction :
		return m_offline.grains().fraction ;
	case sePressure :
		return m_offline.grains().stresses.trace().interpolate< Shape >( m_offline.meshes().primal() ) ;
	case seDivergence :
		return m_offline.grains().sym_grad.trace().interpolate< Shape >( m_offline.meshes().primal() ) ;
	case sePorePressure :
		return m_offline.fluid().pressure ;
	case seFluidVel :
		return m_offline.fluid().velocity.norm().interpolate< Shape >( m_offline.meshes().primal() ) ;
	default:
		return m_offline.grains().spi_grad.interpolate< Shape >( m_offline.meshes().primal() );
	}
}

typename GLViewer::TensorField GLViewer::getTensorEntity() const
{
	switch ( m_tensorEntity ) {
	case teDu:
		return m_offline.grains().sym_grad.interpolate< Shape >( m_offline.meshes().primal() ) ;
	default:
		return m_offline.grains().stresses.interpolate< Shape >( m_offline.meshes().primal() )  ;
	}
}

}
