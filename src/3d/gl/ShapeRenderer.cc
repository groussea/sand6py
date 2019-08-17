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

#include "ShapeRenderer.hh"

#include "MeshRenderer.hh"
#include "TriangularMesh.hh"
#include "Texture.hh"

#include "MeshRenderer.hh"

#include "geo/LevelSet.impl.hh"

namespace d6 {

// Very naive, ugly sphere
static void genSphere( const unsigned parallels, const unsigned meridians,
					   Eigen::Matrix3Xf& ballVerts,
					   std::vector< GLuint >& triIndices )
{
	ballVerts.resize( 3, parallels * meridians ) ;

	triIndices.clear();
	triIndices.reserve( 6 * (parallels - 1) * meridians ) ;

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
				triIndices.push_back( (i-1)*meridians + ( ( j+1 ) % meridians ) ) ;
				triIndices.push_back( (i-1)*meridians + j ) ;
				triIndices.push_back( i*meridians + j ) ;
				triIndices.push_back( (i-1)*meridians + ( ( j+1 ) % meridians ) ) ;
				triIndices.push_back( i*meridians + j ) ;
				triIndices.push_back( i*meridians +  ( ( j+1 ) % meridians ) ) ;
			}
		}
	}

}

static void genPointyCylinder( const float height, const unsigned res,
					   Eigen::Matrix3Xf& vertices,
					   Eigen::Matrix3Xf& normals,
					   Eigen::Matrix3Xf& uvs,
					   std::vector< GLuint >& triIndices )
{
	triIndices.clear() ;
	triIndices.reserve( res * 4 * 3 ) ;
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
	normals.col(1) = Eigen::Vector3f(0,0, 1) ;
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

		triIndices.push_back( 0 ) ;
		triIndices.push_back( id_cur+0 ) ;
		triIndices.push_back( id_nxt+0 ) ;

		triIndices.push_back( id_cur+0 ) ;
		triIndices.push_back( id_nxt+0 ) ;
		triIndices.push_back( id_nxt+1 ) ;
	
		triIndices.push_back( id_cur+0 ) ;
		triIndices.push_back( id_nxt+1 ) ;
		triIndices.push_back( id_cur+1 ) ;

		triIndices.push_back( id_cur+1 ) ;
		triIndices.push_back( id_nxt+1 ) ;
		triIndices.push_back( 1 ) ;

	}

}

static void genTorus( 
				      const float outerRadius,
					  const float radius,
					  const unsigned parallels, const unsigned meridians,
					  Eigen::Matrix3Xf& vertices,
					  Eigen::Matrix3Xf& normals,
					  Eigen::Matrix3Xf& uvs,
					  std::vector< GLuint >& triIndices )
{
	triIndices.clear() ;
	triIndices.reserve( parallels * meridians * 6 ) ;
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

			vertices.col(idx) = outerRadius * p0 + radius * n0 ;

			triIndices.push_back( i*parallels + j );
			triIndices.push_back( i*parallels + ((j+1)%parallels) );
			triIndices.push_back( ((i+1)%meridians)*parallels + ((j+1)%parallels)  );
			
			triIndices.push_back( i*parallels + j );
			triIndices.push_back( ((i+1)%meridians)*parallels + ((j+1)%parallels)  );
			triIndices.push_back( ((i+1)%meridians)*parallels + j  );

		}
	}
}

static void genPlane( const Eigen::Vector3f& box,
					  Eigen::Matrix3Xf& vertices,
					  Eigen::Matrix3Xf& normals,
					  Eigen::Matrix3Xf& uvs,
					  std::vector< GLuint >& triIndices )
{
	triIndices = {0, 1 , 2, 0, 2 , 3} ;
	vertices.resize( 3, 4 ) ;
	normals.resize( 3, 4 ) ;
	uvs.resize( 3, 4 ) ;

	vertices << -box[0], box[0], box[0], -box[0],
		-box[1], -box[1], box[1], box[1],
		0, 0, 0, 0;
	normals.topRows<2>().setZero();
	normals.row(2).setOnes();

	uvs << 0, 1, 1, 0,
		   0, 0, 1, 1,
		   0, 0, 0, 0;
}


static void get_ls_matrix( const LevelSet &ls, Eigen::Matrix4f & mat )
{
	const Eigen::Matrix3f rotation = ls.rotation().matrix().cast < GLfloat >() ;
	const Eigen::Vector3f translation = ls.origin().cast < GLfloat >() ;
	mat.setIdentity() ;
	mat.block<3,3>(0,0) = rotation * ls.scale()  ;
	mat.block<3,1>(0,3) = translation ;
}

static void draw_fake_ball( const LevelSet &ls, const Shader& shader, const gl::ArrayObject& billboardArrays)
{
	const Eigen::Matrix3f rotation = ls.rotation().matrix().cast < GLfloat >() ;
	const Eigen::Vector3f translation = ls.origin().cast < GLfloat >() ;

	glUniform1f( shader.uniform("radius"), ls.scale() ) ;
	glUniformMatrix3fv( shader.uniform("rotation"), 1, GL_FALSE, rotation.data() ) ;
	glUniform3fv( shader.uniform("center"), 1, translation.data() ) ;
	
	//Vertices
	gl::ArrayObject::Using vao(billboardArrays);
	glDrawArrays( GL_TRIANGLES, 0, 6) ;
}

void ShapeRenderer::setup_solid_data(const LevelSet &ls, const Eigen::Vector3f &box,
									 const Shader &shader, ShapeRenderer::MeshDrawData &data)
{
	const CylinderLevelSet *cylinder = dynamic_cast<const CylinderLevelSet *>(&ls);
	const TorusLevelSet *torus = dynamic_cast<const TorusLevelSet *>(&ls);
	const MeshLevelSet *mesh = dynamic_cast<const MeshLevelSet *>(&ls);
	const HoleLevelSet *hole = dynamic_cast<const HoleLevelSet *>(&ls);
	const PlaneLevelSet *plane = dynamic_cast<const PlaneLevelSet *>(&ls);

	Eigen::Matrix3Xf cylVertices, cylNormals, cylUVs;
	std::vector<GLuint> triIndices;

	if (torus)
	{
		genTorus(1.0f, torus->radius(), 30, 30, cylVertices, cylNormals, cylUVs, triIndices);
	}
	else if (cylinder)
	{
		genPointyCylinder(cylinder->height(), 30, cylVertices, cylNormals, cylUVs, triIndices);
	}
	else if (hole)
	{
		genTorus(hole->radius(), 1.0f, 30, 30, cylVertices, cylNormals, cylUVs, triIndices);
	}
	else if (plane)
	{
		genPlane(box, cylVertices, cylNormals, cylUVs, triIndices);
	}

	// Transfer to VBO
	gl::ArrayObject::Using vao(data.vertexArrays);
	data.vertices.reset(cylVertices.cols(), cylVertices.data(), GL_STATIC_DRAW);
	data.normals.reset(cylNormals.cols(), cylNormals.data(), GL_STATIC_DRAW);
	data.uvs.reset(cylUVs.cols(), cylUVs.data(), GL_STATIC_DRAW);
	data.triIndices.reset(triIndices.size(), triIndices.data());

	gl::VertexAttribPointer vap(data.vertices, shader.attribute("vertex"));
	gl::VertexAttribPointer nap(data.normals, shader.attribute("normal"));
	gl::VertexAttribPointer uap(data.uvs, shader.attribute("uv"));
}

void ShapeRenderer::config_shaders()
{
	// Ball shader
	m_ballShader.add_uniform("model_view") ;
	m_ballShader.add_uniform("projection") ;
	m_ballShader.add_uniform("radius") ;
	m_ballShader.add_uniform("rotation") ;
	m_ballShader.add_uniform("center") ;
	m_ballShader.add_uniform("light_pos") ;
	m_ballShader.add_attribute("vertex") ;
	m_ballShader.load("ball_vertex","ball_fragment") ;
	
	// Ball depth shader
	m_ballDepthShader.add_uniform("model_view") ;
	m_ballDepthShader.add_uniform("projection") ;
	m_ballDepthShader.add_uniform("radius") ;
	m_ballDepthShader.add_uniform("rotation") ;
	m_ballDepthShader.add_uniform("center") ;
	m_ballDepthShader.add_attribute("vertex") ;
	m_ballDepthShader.load("ball_depth_vertex","ball_depth_fragment") ;

	// Default solid shader
	m_solidShader.add_uniform("model_view") ;
	m_solidShader.add_uniform("projection") ;
	m_solidShader.add_uniform("light_pos") ;
	m_solidShader.add_uniform("ambient") ;
	m_solidShader.add_uniform("depth_mvp") ;
	m_solidShader.add_uniform("depth_texture") ;
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

void ShapeRenderer::setup_vaos()
{

	// Gen glyph vertices
	{
		Eigen::Matrix3Xf sphereVertices;
		std::vector<GLuint> triIndices;
		genSphere(5, 8, sphereVertices, triIndices);
		m_sphereVertices.reset(sphereVertices.cols(), sphereVertices.data(), GL_STATIC_DRAW);
		
		// Bind tri indices to array object
		gl::ArrayObject::Using vao(m_sphereVertexArrays);
		m_sphereTriIndices.reset( triIndices.size(), triIndices.data() );
	}

	{
		Eigen::Matrix<float, 3, 6> vtx;
		vtx << -1, -1, 1, -1, 1, 1,
			-1, 1, 1, -1, 1, -1,
			0, 0, 0, 0, 0, 0;
		gl::ArrayObject::Using vao(m_billboardArrays);
		m_squareVertices.reset(vtx.cols(), vtx.data());

		gl::VertexAttribPointer vap_v(m_squareVertices, m_ballShader.attribute("vertex"));
	}

}

void ShapeRenderer::init()
{
	config_shaders();
	setup_vaos();
}

void ShapeRenderer::clear_buffers()
{
	m_solidData.clear();
}
	
void ShapeRenderer::setup_buffers(const LevelSet &ls, const Eigen::Vector3f& box)
{
	if( dynamic_cast<const SphereLevelSet*>(&ls) ) return;

	setup_solid_data(ls, box, m_solidShader, m_solidData[&ls]);
}

void ShapeRenderer::compute_shadow( const LevelSet &ls, 
		const Eigen::Matrix4f& depthModelView, const Eigen::Matrix4f& depthProjection ) const 
{

	if( dynamic_cast<const SphereLevelSet*>(&ls) )
	{
		UsingShader sh( m_ballDepthShader ) ;
		sh.bindMVP(depthModelView.data(), depthProjection.data());
	
		draw_fake_ball( ls, m_ballDepthShader, m_billboardArrays ) ;
	} else {
		// Solid shader
		Eigen::Matrix4f mat ;
		get_ls_matrix( ls, mat ) ;

		Eigen::Matrix4f completeMVP = depthProjection * depthModelView * mat ;

		UsingShader sh( m_solidDepthShader ) ;
		glUniformMatrix4fv( m_solidDepthShader.uniform("depth_mvp"), 1, GL_FALSE, completeMVP.data()) ;

		auto dataIt = m_solidData.find(&ls);
		if (dataIt != m_solidData.end())
		{
			const auto &data = dataIt->second;
			gl::ArrayObject::Using vao(data.vertexArrays);
			glDrawElements(GL_TRIANGLES, data.triIndices.size(), GL_UNSIGNED_INT, 0);
		}
	}
}

void ShapeRenderer::draw( const LevelSet &ls, const Vec &box, const Eigen::Vector3f& lightPos,
		bool shadowed, const Texture& depthTexture, 
		const Eigen::Matrix4f& modelView, const Eigen::Matrix4f& projection,
		const Eigen::Matrix4f& depthModelView, const Eigen::Matrix4f& depthProjection ) const 
{
	const Eigen::Matrix3f rotation = ls.rotation().matrix().cast < GLfloat >() ;
	const Eigen::Vector3f translation = ls.origin().cast < GLfloat >() ;

	if( dynamic_cast<const SphereLevelSet*>(&ls) )
	{
		UsingShader sh( m_ballShader ) ;
		// Model-view
		sh.bindMVP(modelView.data(), projection.data()) ;

		draw_fake_ball( ls, m_ballShader, m_billboardArrays ) ;
	} else {
		Eigen::Matrix4f mat ;
		get_ls_matrix( ls, mat ) ;
		Eigen::Matrix4f completeDepthMVP = depthProjection * depthModelView * mat ;

		//glColor4f(1., 0., .8, 1);
		Eigen::Matrix4f finalModelView = modelView * mat; 

		const Eigen::Vector3f objLight = rotation.inverse() / ls.scale() * (lightPos - translation) ;

			auto dataIt = m_solidData.find(&ls);
			Eigen::Vector3f color(0.3, 0.15, 0.1);

			if (m_solidShader.ok() && dataIt != m_solidData.end())
			{
				UsingShader sh(m_solidShader);
				sh.bindMVP(finalModelView.data(), projection.data());

				glUniform3fv(m_solidShader.uniform("light_pos"), 1, objLight.data());
				glUniform3fv(m_solidShader.uniform("ambient"), 1, color.data());

				UsingTexture tx(depthTexture);
				if (shadowed)
				{
					tx.bindUniform(m_solidShader.uniform("depth_texture"));
					glUniformMatrix4fv(m_solidShader.uniform("depth_mvp"), 1, GL_FALSE, completeDepthMVP.data());
				}

				const auto &data = dataIt->second;
				gl::ArrayObject::Using vao(data.vertexArrays);
				glDrawElements(GL_TRIANGLES, data.triIndices.size(), GL_UNSIGNED_INT, 0);

			}
	}


}

} // d6
