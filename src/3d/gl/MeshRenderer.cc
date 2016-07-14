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

#include "MeshRenderer.hh"
#include "Shader.hh"

#include "geo/TriangularMesh.hh"

#include <iostream>

namespace d6 {

void MeshRenderer::reset( const TriangularMesh& mesh )
{

	const unsigned n = mesh.nFaces() ;

	Eigen::Matrix3Xf vertices( 3, n*3 ) ;

#pragma omp parallel for
	for( unsigned i = 0 ; i < n ; ++i ) {
		vertices.col(3*i+0) = mesh.vertex( i,0 ).cast<float>() ;
		vertices.col(3*i+1) = mesh.vertex( i,1 ).cast<float>() ;
		vertices.col(3*i+2) = mesh.vertex( i,2 ).cast<float>() ;
	}
	m_vertices.reset( 3*n, vertices.data() ) ;

	if ( mesh.hasVertexNormals() ) {
		Eigen::Matrix3Xf normals( 3, n*3 ) ;

#pragma omp parallel for
		for( unsigned i = 0 ; i < n ; ++i ) {
			normals.col(3*i+0) = mesh.normal( i,0 ).cast<float>() ;
			normals.col(3*i+1) = mesh.normal( i,1 ).cast<float>() ;
			normals.col(3*i+2) = mesh.normal( i,2 ).cast<float>() ;
		}
		m_normals.reset( 3*n, normals.data() ) ;

	} else if ( mesh.hasFaceNormals() ) {
		Eigen::Matrix3Xf normals( 3, n*3 ) ;

#pragma omp parallel for
		for( unsigned i = 0 ; i < n ; ++i ) {
			normals.col(3*i+0) = mesh.faceNormal( i ).cast<float>() ;
			normals.col(3*i+1) = mesh.faceNormal( i ).cast<float>() ;
			normals.col(3*i+2) = mesh.faceNormal( i ).cast<float>() ;
		}
		m_normals.reset( 3*n, normals.data() ) ;

	}

	if( mesh.hasVertexUVs() )
	{
		Eigen::Matrix3Xf uvs( 3, n*3 ) ;

#pragma omp parallel for
		for( unsigned i = 0 ; i < n ; ++i ) {
			uvs.col(3*i+0) = mesh.uv( i,0 ).cast<float>() ;
			uvs.col(3*i+1) = mesh.uv( i,1 ).cast<float>() ;
			uvs.col(3*i+2) = mesh.uv( i,2 ).cast<float>() ;
		}
		m_uvs.reset( 3*n, uvs.data() ) ;
	}

}

void MeshRenderer::draw( const Shader &shader ) const
{
	gl::VertexPointer vp( m_vertices ) ;
	if( shader.ok() )
	{
		gl::VertexAttribPointer vap( m_vertices, shader.attribute("vertex") ) ;
		gl::VertexAttribPointer nap( m_normals,  shader.attribute("normal") ) ;
		gl::VertexAttribPointer uap( m_uvs,  shader.attribute("uv") ) ;

		glDrawArrays( GL_TRIANGLES, 0, m_vertices.size() ) ;

	} else {
		gl::NormalPointer np( m_normals ) ;

		glDrawArrays( GL_TRIANGLES, 0, m_vertices.size() ) ;
	}
}

} //d6

