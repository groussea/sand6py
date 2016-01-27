#include "GLViewer.hh"

#include "geo/MeshImpl.hh"

#include "simu/Phase.hh"
#include "visu/Offline.hh"

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


void GLViewer::init( )
{
	updateViewport();


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


void GLViewer::update_buffers()
{
	const MeshType& g = m_offline.mesh() ;

	// Grid nodes
	// Grid quad indices
	{
		Eigen::Matrix<float, 2, Eigen::Dynamic > vertices( 2, g.nNodes() ) ;
		std::vector<unsigned> nodeIndices( MeshType::NV * g.nCells() ) ;

		typename MeshType::CellGeo cellGeo ;
		typename MeshType::NodeList cellNodes ;

		for( typename MeshType::CellIterator it = g.cellBegin() ; it != g.cellEnd() ; ++it )
		{
			g.get_geo( *it, cellGeo ) ;
			g.list_nodes( *it, cellNodes );

			nodeIndices[ 4*it.index() + 0 ] = cellNodes[0] ;
			nodeIndices[ 4*it.index() + 1 ] = cellNodes[1] ;
			nodeIndices[ 4*it.index() + 2 ] = cellNodes[3] ;
			nodeIndices[ 4*it.index() + 3 ] = cellNodes[2] ;

			for( int k = 0 ; k < MeshType::NV ; ++k ) {
				vertices.col( cellNodes[k] ) = cellGeo.vertex( k ).cast< float >() ;
			}
		}

		m_gridVertices.reset( vertices.cols(), vertices.data() ) ;
		m_gridQuadIndices.reset( nodeIndices.size(), &nodeIndices[0] ) ;
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
	std::vector< float > densities ( p.count() );
	std::vector< float > velocities ( p.count() );
	Eigen::Matrix< float, 4, Eigen::Dynamic > frames (4, p.count()) ;

	// Compute movel-view matrix from tensor
	Mat tensor ;
	Eigen::Matrix2f mat ;

#pragma omp parallel for private(tensor)
	for( size_t i = 0 ; i < p.count() ; ++i ) {

		if( m_drawOrientations) {
			tensor_view( p.orient().col( i ) ).get( tensor ) ;
			Eigen::SelfAdjointEigenSolver<Mat> es( tensor );

			const Vec ev = es.eigenvalues().array().max(0).sqrt() ;

			mat = ( es.eigenvectors() * ev.asDiagonal() ).cast< GLfloat >()
					* .5 * std::pow( p.volumes()[i], 1./3 )  ;

			densities[i] = p.volumes()[i]  ;

		} else {
			tensor_view( p.frames().col( i ) ).get( tensor ) ;
			Eigen::SelfAdjointEigenSolver<Mat> es( tensor );

			const Vec ev = es.eigenvalues().array().max(0).sqrt() ;
			const Scalar vol = 8 * ev.prod() ;

			mat = ( es.eigenvectors() * ev.asDiagonal() ).cast< GLfloat >()  ;

			densities[i] = p.volumes()[i] / vol ;
		}

		velocities[i] = p.centers().col(i).norm()  ;


		frames.col(i) = Eigen::Matrix< GLfloat, 4, 1 >::Map( mat.data(), mat.size() ) ;

	}
	m_particleFrames.reset( p.count(), frames.data(), GL_STATIC_DRAW )  ;
	m_particleAlpha.reset ( p.count(), densities.data(), GL_STATIC_DRAW )  ;
	m_particleColors.reset ( p.count(), velocities.data(), GL_STATIC_DRAW )  ;
}

// Grid colors
void GLViewer::update_color_buffers()
{
	const MeshType& g = m_offline.mesh() ;

	Eigen::Matrix<float, 3, Eigen::Dynamic > colors( 3, g.nNodes() ) ;
	colors.setZero() ;

	const auto & field = getScalarEntity() ;

	const double max_val = field.flatten().maxCoeff() ;
	const double min_val = field.flatten().minCoeff() ;
	if( max_val - min_val > 1.e-12 )
	{
		colors.row(1) = ( (field.flatten().array() - min_val)/(max_val - min_val) ).cast<float>() ;
	}

	m_gridColors.reset( colors.cols(), colors.data() ) ;
}


// Arrows
void GLViewer::update_vector_buffers()
{
	const MeshType& g = m_offline.mesh() ;

	Eigen::Matrix<float, 2, Eigen::Dynamic > vectors( 2, g.nNodes() ) ;
	Eigen::Matrix<float, 3, Eigen::Dynamic >  colors( 3, g.nNodes() ) ;

	vectors.setZero() ;
	colors.setZero() ;

	const auto & field = getVectorEntity() ;
	ScalarField norm = field.norm() ;

	const double max_val = norm.max_abs() ;
	if( max_val > 1.e-12 )
	{
		Eigen::VectorXf::Map( vectors.data(), vectors.size() ) = ( field.flatten().array() / max_val ).cast<float>() ;
		colors.row(2) =  ( norm.flatten().array() / max_val ).cast<float>() ;
	}
	colors.row(1).setOnes() ;

	m_arrows.reset( vectors.cols(), vectors.data() ) ;
	m_arrowColors.reset( colors.cols(), colors.data() ) ;
}

// Arrows
void GLViewer::update_tensor_buffers()
{
	const MeshType& g = m_offline.mesh() ;

	Eigen::Matrix<float, 4, Eigen::Dynamic > tensors( 2, g.nNodes() ) ;
	Eigen::Matrix<float, 3, Eigen::Dynamic >  colors( 3, g.nNodes() ) ;

	tensors.setZero() ;
	colors.setZero() ;

	const auto & field = getTensorEntity() ;
	ScalarField norm = field.norm() ;

	const double max_val = norm.max_abs() ;
	if( max_val > 1.e-12 )
	{
		Eigen::VectorXf::Map( tensors.data(), tensors.size() ) = ( field.flatten().array() / max_val ).cast<float>() ;
		colors.row(2) =  ( norm.flatten().array() / max_val ).cast<float>() ;
	}
	colors.row(1).setOnes() ;
	colors.row(0).array() = 1 -  colors.row(2).array();

	m_tensors.reset( tensors.cols(), tensors.data() ) ;
	m_tensorColors.reset( colors.cols(), colors.data() ) ;
}

void GLViewer::update_texture()
{
	const Config& c = m_offline.config() ;

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

			Vec u = field(Vec(x*sx,y*sy)) ;
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
				u = field(Vec(x*sx,y*sy)) ;
			}
			x = c ; y = r ;
			for( unsigned k = 0 ; k < kmax ; ++k )
			{
				u = field(Vec(x*sx,y*sy)) ;
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

void GLViewer::draw( )
{
	const Config& c = m_offline.config() ;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glLineWidth( 1.f );

	// Grid

	if( m_shouldRender[ eColors ])
	{
		m_gridQuadIndices.bind() ;
		gl::VertexPointer vp( m_gridVertices ) ;
		gl::ColorPointer cp( m_gridColors ) ;
		glDrawElements( GL_QUADS, m_gridQuadIndices.size(), GL_UNSIGNED_INT, 0 );
	}

	glColor3f( .5, .5, .5 );
	if( m_shouldRender[ eGrid ])
	{
		m_gridQuadIndices.bind() ;
		gl::VertexPointer vp( m_gridVertices ) ;
		glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
		glDrawElements( GL_QUADS, m_gridQuadIndices.size(), GL_UNSIGNED_INT, 0 );
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

		glColor3f( 1., 0., 0. );
		glPointSize( 2.f ) ;

		gl::VertexPointer vp( m_particles ) ;
		glDrawArrays( GL_POINTS, 0, m_particles.size() );
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


//	//LS
//	for( unsigned i = 0 ; i < m_simu.nLevelSets() ; ++i ) {
//		glColor3f( 0., 0., 1. );
//		glLineWidth( 2.f ) ;
//		m_simu.levelSet(i).draw();
//	}
}

const VectorField& GLViewer::getVectorEntity() const
{
	switch ( m_vectorEntity ) {
	case veVelocity :
		return m_offline.grains().velocity ;
	case veForce :
		return m_offline.grains().fcontact ;
	case veProj :
		return m_offline.grains().geo_proj ;
	default:
		return m_offline.grains().grad_phi ;
	}
}
ScalarField GLViewer::getScalarEntity() const
{
	switch ( m_scalarEntity ) {
	case seFraction :
		return m_offline.grains().fraction ;
	case sePressure :
		return m_offline.grains().stresses.trace() ;
	default:
		return m_offline.grains().spi_grad ;
	}
}

const TensorField& GLViewer::getTensorEntity() const
{
	switch ( m_tensorEntity ) {
	case teDu:
		return m_offline.grains().sym_grad ;
	default:
		return m_offline.grains().stresses ;
	}
}

}
