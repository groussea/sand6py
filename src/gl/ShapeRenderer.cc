#include "ShapeRenderer.hh"

#include "geo/LevelSet.impl.hh"

namespace d6 {

// Very naive, ugly sphere
static void genSphere( const unsigned parallels, const unsigned meridians,
					   Eigen::Matrix3Xf& ballVerts,
					   std::vector< GLuint >& quadIndices )
{
	ballVerts.resize( 3, parallels * meridians ) ;

	quadIndices.clear();
	quadIndices.reserve( 4 * (parallels - 1) * meridians ) ;

	const double dphi = M_PI / (parallels - 1.);
	const double dpsi = 2 * M_PI / meridians ;

	for(unsigned i = 0; i < parallels; ++i)
	{
		const double z = std::cos( i * dphi ) ;
		const double r = std::sin( i * dphi ) ;

		for(unsigned j = 0; j < meridians; j++)
		{
			const Vec p ( r * std::cos( j * dpsi ), r * std::sin( j * dpsi ), z )  ;
			ballVerts.col( i*meridians + j ) = p.cast< GLfloat >() ;

			if( i > 0 ) {
				quadIndices.push_back( (i-1)*meridians + ( ( j+1 ) % meridians ) ) ;
				quadIndices.push_back( (i-1)*meridians + j ) ;
				quadIndices.push_back( i*meridians + j ) ;
				quadIndices.push_back( i*meridians +  ( ( j+1 ) % meridians ) ) ;
			}
		}
	}

}

void ShapeRenderer::init()
{
	// Gen glyph vertices
	Eigen::Matrix3Xf sphereVertices ;
	std::vector< GLuint > quadIndices ;
	genSphere( 5, 8, sphereVertices, quadIndices );
	m_sphereVertices.reset( sphereVertices.cols(), sphereVertices.data(), GL_STATIC_DRAW );
	m_sphereQuadIndices.reset( quadIndices.size(), quadIndices.data() );

	Eigen::Matrix<float, 3, 4> vtx ;
	vtx  <<  -1, -1, 1,  1,
		 -1,  1, 1, -1,
		 0, 0 ,0, 0 ;
	m_squareVertices.reset( 4, vtx.data() );


	// Ball shader
	m_ballShader.add_uniform("model_view") ;
	m_ballShader.add_uniform("projection") ;
	m_ballShader.add_uniform("radius") ;
	m_ballShader.add_uniform("rotation") ;
	m_ballShader.add_uniform("center") ;
	m_ballShader.add_uniform("light_pos") ;
	m_ballShader.add_attribute("vertex") ;
	m_ballShader.load("ball_vertex","ball_fragment") ;

}

void ShapeRenderer::draw( const LevelSet &ls, const Vec &box, const Eigen::Vector3f& lightPos ) const
{
	const Eigen::Matrix3f rotation = ls.rotation().matrix().cast < GLfloat >() ;
	const Eigen::Vector3f translation = ls.origin().cast < GLfloat >() ;

	if( dynamic_cast<const SphereLevelSet*>(&ls) )
	{

		UsingShader sh( m_ballShader ) ;
		// Model-view
		sh.bindMVP() ;

		//Vertices
		gl::VertexAttribPointer vap_v( m_squareVertices, m_ballShader.attribute("vertex") ) ;

		glUniform1f( m_ballShader.uniform("radius"), ls.scale() ) ;
		glUniformMatrix3fv( m_ballShader.uniform("rotation"), 1, GL_FALSE, rotation.data() ) ;
		glUniform3fv( m_ballShader.uniform("center"), 1, translation.data() ) ;

		glUniform3fv( m_ballShader.uniform("light_pos"), 1, lightPos.data() ) ;

		gl::VertexPointer vp( m_squareVertices ) ;
		glDrawArrays( GL_QUADS, 0, m_squareVertices.size() ) ;

	} else {

		Eigen::Matrix4f mat ;
		mat.setIdentity() ;
		mat.block<3,3>(0,0) = rotation * ls.scale()  ;
		mat.block<3,1>(0,3) = translation ;

		glColor4f(1., .8, .8, 1);


		glPushMatrix();
		glMultMatrixf( mat.data() );

		const CylinderLevelSet* cylinder = nullptr ;
		const TorusLevelSet* torus = nullptr ;

		if ( dynamic_cast<const PlaneLevelSet*>(&ls) ) {

			glBegin( GL_QUADS );
			glNormal3f( 0.f, 0.f, 1.f );
			glVertex3d( -box[0], -box[1], 0 );
			glNormal3f( 0.f, 0.f, 1.f );
			glVertex3d( -box[0],  box[1], 0 );
			glNormal3f( 0.f, 0.f, 1.f );
			glVertex3d(  box[0],  box[1], 0 );
			glNormal3f( 0.f, 0.f, 1.f );
			glVertex3d(  box[0], -box[1], 0 );
			glEnd( ) ;

		} else if ( (torus = dynamic_cast<const TorusLevelSet*>(&ls)) ) {

			const unsigned res= 10 ;
			for( unsigned i = 0 ; i < res ; ++i ) {

				const float alpha0 = (2*M_PI*(i+0)) / res ;
				const float alpha1 = (2*M_PI*(i+1)) / res ;

				const Eigen::Vector3f p0( std::cos(alpha0), std::sin(alpha0), 0 ) ;
				const Eigen::Vector3f p1( std::cos(alpha1), std::sin(alpha1), 0 ) ;

				glBegin( GL_QUAD_STRIP );
				for( unsigned j = 0 ; j < res ; ++j )
				{
					const float beta0 = (2*M_PI*(j+0)) / (res-1) ;

					Eigen::Vector3f n0 = std::cos(beta0)*p0.normalized() ;
					n0[2] = std::sin(beta0) ;

					Eigen::Vector3f n1 = std::cos(beta0)*p1.normalized() ;
					n1[2] = std::sin(beta0) ;

					Eigen::Vector3f v0 = p0 + torus->radius() * n0 ;
					Eigen::Vector3f v1 = p1 + torus->radius() * n1 ;

					glNormal3fv( n0.data() );
					glVertex3fv( v0.data() );
					glNormal3fv( n1.data() );
					glVertex3fv( v1.data() );

				}
				glEnd( ) ;
			}

		} else if ( (cylinder = dynamic_cast<const CylinderLevelSet*>(&ls)) ) {

		}

		glPopMatrix();
	}

}

} // d6
