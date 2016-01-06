#include "ShapeRenderer.hh"

#include "MeshRenderer.hh"
#include "TriangularMesh.hh"

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

static void genPointyCylinder( const float height, const unsigned res,
					   Eigen::Matrix3Xf& vertices,
					   Eigen::Matrix3Xf& normals,
					   Eigen::Matrix3Xf& uvs,
					   std::vector< GLuint >& quadIndices )
{
	quadIndices.clear() ;
	quadIndices.reserve( res * 4 * 3 ) ;
	vertices.resize( 3, 2 + res * 2 ) ;
	normals.resize( 3, 2 + res * 2 ) ;
	uvs.resize( 3, 2 + res * 2 ) ;

	const double dphi = 2*M_PI / (res);

	Eigen::Vector3f p0 (0,0,-.5*height) ;
	Eigen::Vector3f p1 (0,0, .5*height) ;

	vertices.col(0) = p0 + Eigen::Vector3f(0,0,-1) ;
	normals.col(0) = Eigen::Vector3f(0,0,-1) ;
	uvs.col(0) = Eigen::Vector3f(0,0,0) ;
	vertices.col(1) = p1 + Eigen::Vector3f(0,0, 1);
	uvs.col(1) = Eigen::Vector3f(1,0, 0) ;

	for(unsigned i = 0; i < res; ++i)
	{
		const float beta0 = i*dphi ;
		const unsigned id_cur = 2+2*i ;
		const unsigned id_nxt = 2+2*((i+1)%res) ;

		Eigen::Vector3f n ( std::cos(beta0), std::sin(beta0),0) ;
		vertices.col(id_cur+0) = p0 + n ;
		vertices.col(id_cur+1) = p1 + n ;
		normals .col(id_cur+0) = n ;
		normals .col(id_cur+1) = n ;
		uvs     .col(id_cur+0) = Eigen::Vector3f(1/(height+2),beta0/(2*M_PI),0) ;
		uvs     .col(id_cur+1) = Eigen::Vector3f((1+height)/(height+2),beta0/(2*M_PI),0) ;

		quadIndices.push_back( 0 ) ;
		quadIndices.push_back( id_cur+0 ) ;
		quadIndices.push_back( id_nxt+0 ) ;
		quadIndices.push_back( 0 ) ;

		quadIndices.push_back( id_cur+0 ) ;
		quadIndices.push_back( id_nxt+0 ) ;
		quadIndices.push_back( id_nxt+1 ) ;
		quadIndices.push_back( id_cur+1 ) ;

		quadIndices.push_back( 1 ) ;
		quadIndices.push_back( id_cur+1 ) ;
		quadIndices.push_back( id_nxt+1 ) ;
		quadIndices.push_back( 1 ) ;

	}

}

static void genTorus( const float radius,
					  const unsigned parallels, const unsigned meridians,
					  Eigen::Matrix3Xf& vertices,
					  Eigen::Matrix3Xf& normals,
					  Eigen::Matrix3Xf& uvs,
					  std::vector< GLuint >& quadIndices )
{
	quadIndices.clear() ;
	quadIndices.reserve( parallels * meridians * 4 ) ;
	vertices.resize( 3, parallels * meridians ) ;
	normals.resize( 3, parallels * meridians ) ;
	uvs.resize( 3, parallels * meridians ) ;

	const Scalar dalpha = (2*M_PI)/(meridians-1) ;
	const Scalar dbeta  = (2*M_PI)/(parallels-1) ;

	for( unsigned i = 0 ; i < meridians ; ++i )
	{

		const float alpha0 = i*dalpha  ;
		const Eigen::Vector3f p0 ( std::cos(alpha0), std::sin(alpha0), 0 ) ;

		for( unsigned j = 0 ; j < parallels ; ++j )
		{
			const unsigned idx = i*parallels + j ;
			const float beta0 = j * dbeta ;

			Eigen::Vector3f n0 = std::cos(beta0)*p0.normalized() ;
			n0[2] = std::sin(beta0) ;
			normals.col(idx) = n0 ;
			uvs.col(idx) = Eigen::Vector3f(alpha0/(2*M_PI),beta0/(2*M_PI),0) ;

			vertices.col(idx) = p0 + radius * n0 ;

			quadIndices.push_back( i*parallels + j );
			quadIndices.push_back( i*parallels + ((j+1)%parallels) );
			quadIndices.push_back( ((i+1)%meridians)*parallels + ((j+1)%parallels)  );
			quadIndices.push_back( ((i+1)%meridians)*parallels + j  );

		}
	}
}

static void get_ls_matrix( const LevelSet &ls, Eigen::Matrix4f & mat )
{
	const Eigen::Matrix3f rotation = ls.rotation().matrix().cast < GLfloat >() ;
	const Eigen::Vector3f translation = ls.origin().cast < GLfloat >() ;
	mat.setIdentity() ;
	mat.block<3,3>(0,0) = rotation * ls.scale()  ;
	mat.block<3,1>(0,3) = translation ;
}

static void draw_solid( const Shader& shader, 
			const Eigen::Matrix3Xf& cylVertices,
			const Eigen::Matrix3Xf& cylNormals,
			const Eigen::Matrix3Xf& cylUVs,
			const std::vector< GLuint > &quadIndices )
{
	// Transfer to VBO
	gl::VertexBuffer3f cylv, cyln, cyluv ;
	gl::IndexBuffer cylqi ;
	cylv.reset( cylVertices.cols(), cylVertices.data(), GL_DYNAMIC_DRAW );
	cyln.reset( cylNormals .cols(), cylNormals .data(), GL_DYNAMIC_DRAW );
	cyluv.reset( cylUVs    .cols(), cylUVs     .data(), GL_DYNAMIC_DRAW );
	cylqi.reset( quadIndices.size(), quadIndices.data() );

	cylqi.bind() ;
	
	gl::VertexAttribPointer vap( cylv, shader.attribute("vertex") ) ;
	gl::VertexAttribPointer nap( cyln, shader.attribute("normal") ) ;
	gl::VertexAttribPointer uap( cyln, shader.attribute("uv") ) ;

	glDrawElements( GL_QUADS, cylqi.size(), GL_UNSIGNED_INT, 0 );
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

	// Default solid shader
	m_solidShader.add_uniform("model_view") ;
	m_solidShader.add_uniform("projection") ;
	m_solidShader.add_uniform("light_pos") ;
	m_solidShader.add_uniform("ambient") ;
	m_solidShader.add_attribute("vertex") ;
	m_solidShader.add_attribute("normal") ;
	m_solidShader.add_attribute("uv") ;
	m_solidShader.load("vertex","fragment") ;
	
	// Default solid depth shader
	m_solidDepthShader.add_uniform("depth_mvp") ;
	m_solidDepthShader.add_attribute("vertex") ;
	m_solidDepthShader.add_attribute("normal") ;
	m_solidDepthShader.add_attribute("uv") ;
	m_solidDepthShader.load("depth_vertex","depth_fragment") ;

}

void ShapeRenderer::compute_shadow( const LevelSet &ls, const Eigen::Matrix4f& depthMVP ) const
{

	if( dynamic_cast<const SphereLevelSet*>(&ls) )
	{
		// TODO
	} else {
		// Solid shader
		Eigen::Matrix4f mat ;
		get_ls_matrix( ls, mat ) ;
	
		Eigen::Matrix4f completeMVP = depthMVP * mat ;

		const CylinderLevelSet* cylinder = nullptr ;
		const TorusLevelSet* torus = nullptr ;
		//const HoleLevelSet* hole = nullptr ;
		//const MeshLevelSet* mesh = nullptr ;

		if ( (torus = dynamic_cast<const TorusLevelSet*>(&ls))
				 || (cylinder = dynamic_cast<const CylinderLevelSet*>(&ls)) ) {

			Eigen::Matrix3Xf cylVertices, cylNormals, cylUVs ;
			std::vector< GLuint > quadIndices ;

			if( torus ) {
				genTorus( torus->radius(), 30, 30, cylVertices, cylNormals, cylUVs, quadIndices );
			} else {
				genPointyCylinder( cylinder->height(), 30, cylVertices, cylNormals, cylUVs, quadIndices );
			}

			UsingShader sh( m_solidDepthShader ) ;
			glUniformMatrix4fv( m_solidDepthShader.uniform("depth_mvp"), 1, GL_FALSE, completeMVP.data()) ;

			draw_solid( m_solidDepthShader, cylVertices, cylNormals, cylUVs, quadIndices ) ;
		}	
	}
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
		get_ls_matrix( ls, mat ) ; 

		glColor4f(1., 0., .8, 1);

		glPushMatrix();
		glMultMatrixf( mat.data() );

		const Eigen::Vector3f objLight = rotation.inverse() / ls.scale() * (lightPos - translation) ;

		const CylinderLevelSet* cylinder = nullptr ;
		const TorusLevelSet* torus = nullptr ;
		const HoleLevelSet* hole = nullptr ;
		const MeshLevelSet* mesh = nullptr ;

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

		} else if ( (hole = dynamic_cast<const HoleLevelSet*>(&ls)) ) {

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

					Eigen::Vector3f v0 = hole->radius() * p0 + n0 ;
					Eigen::Vector3f v1 = hole->radius() * p1 + n1 ;

					glNormal3fv( n0.data() );
					glVertex3fv( v0.data() );
					glNormal3fv( n1.data() );
					glVertex3fv( v1.data() );

				}
				glEnd( ) ;
			}
		} else if ( (torus = dynamic_cast<const TorusLevelSet*>(&ls))
				 || (cylinder = dynamic_cast<const CylinderLevelSet*>(&ls)) ) {

			// Generate points = normales
			Eigen::Matrix3Xf cylVertices, cylNormals, cylUVs ;
			std::vector< GLuint > quadIndices ;

			if( torus ) {
				genTorus( torus->radius(), 30, 30, cylVertices, cylNormals, cylUVs, quadIndices );
			} else {
				genPointyCylinder( cylinder->height(), 30, cylVertices, cylNormals, cylUVs, quadIndices );
			}

			Eigen::Vector3f color(0.3,0.15,0.1) ;

			if( m_solidShader.ok() ) {
				UsingShader sh( m_solidShader ) ;
				sh.bindMVP() ;

				glUniform3fv( m_solidShader.uniform("light_pos"), 1, objLight.data() ) ;
				glUniform3fv( m_solidShader.uniform("ambient"), 1, color.data() ) ;

				draw_solid( m_solidShader, cylVertices, cylNormals, cylUVs, quadIndices ) ;
			}
		} else if ( (mesh = dynamic_cast<const MeshLevelSet*>(&ls)) ) {
			MeshRenderer& renderer = m_meshRenderers[ mesh->objFile() ] ;
			if( !renderer.ok() ) {
				TriangularMesh triMesh ;
				triMesh.loadObj( mesh->objFile().c_str() ) ;
				if( !triMesh.hasVertexNormals() )
					triMesh.computeFaceNormals() ;

				renderer.reset( triMesh ) ;
			}

			UsingShader sh( m_solidShader ) ;

			Eigen::Vector3f color(0.3,0.15,0.1) ;
			glColor3fv( color.data() ) ;

			if( m_solidShader.ok() ) {
				sh.bindMVP() ;
				glUniform3fv( m_solidShader.uniform("light_pos"), 1, objLight.data() ) ;
				glUniform3fv( m_solidShader.uniform("ambient"), 1, color.data() ) ;
			}

			renderer.draw( m_solidShader ) ;
		}


		glPopMatrix();
	}


}

} // d6
