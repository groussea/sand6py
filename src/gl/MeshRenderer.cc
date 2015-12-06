#include "MeshRenderer.hh"

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

}

void MeshRenderer::draw() const {
	
	gl::VertexPointer vp( m_vertices ) ;
	gl::NormalPointer np( m_normals ) ;

	glDrawArrays( GL_TRIANGLES, 0, m_vertices.size() ) ;

}

} //d6

