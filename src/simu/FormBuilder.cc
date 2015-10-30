#include "FormBuilder.hh"

#include "geo/Grid.hh"

#include <bogus/Core/Block.impl.hpp>
#include <algorithm>

namespace d6 {

void FormBuilder::addToIndex(
		const typename MeshType::Cells& cells,
		const std::vector< Index > &rowIndices,
		const std::vector< Index > &colIndices
		) {

	typename MeshType::NodeList nodes ;
	for( const typename MeshType::Cell& cell : cells ) {
		m_mesh.list_nodes( cell, nodes );
		for( int k = 0 ; k < MeshType::NV ; ++ k ) {
			for( int j = 0 ; j < MeshType::NV ; ++ j ) {
				m_data[ rowIndices[ nodes[k] ] ].push_back( colIndices[ nodes[j] ] ) ;
			}
		}
	}

}

void FormBuilder::makeCompressed()
{

	const size_t m = m_data.size() ;

#pragma omp parallel for
	for( size_t i = 0 ; i < m ; ++i ) {
		std::sort( m_data[i].begin(),  m_data[i].end() ) ;
		auto last = std::unique( m_data[i].begin(),  m_data[i].end() )  ;
		m_data[i].erase( last,  m_data[i].end() ) ;
	}

	m_compressed.clear() ;
	m_compressed.outer.resize( m + 1 ) ;

	m_compressed.outer[0] = 0 ;
	for( size_t i = 0 ; i < m ; ++i ) {
		m_compressed.outer[i+1] = m_compressed.outer[i] + m_data[i].size() ;
	}

	m_compressed.inner.resize( m_compressed.outer.back() ) ;

#pragma omp parallel for
	for( size_t i = 0 ; i < m ; ++i ) {
		memcpy( m_compressed.inner.data() + m_compressed.outer[i], m_data[i].data(),
				m_compressed.size( i ) * sizeof( BgIndex ) ) ;
	}
}

/*
 *
 * D( u )
 * a = Tr D(u) / sqrt(6) = ( dux_dx + duy_dy + duz_dz ) / sqrt(6)
 * b = .5 ( dux_dx - duy_dy )
 * c = ( 2*duz_dz - duy_dy - dux_dx ) / (2 * s_sqrt_3)
 * d = .5 (dux_dy + duy_dx)
 * e = .5 (dux_dz + duz_dx)
 * f = .5 (duz_dy + duy_dz)
 */

void FormBuilder::addDuDv( FormMat<3,3>::Type& A, Scalar w, Itp itp, Dcdx dc_dx, Indices rowIndices, Indices colIndices )
{
	typedef FormMat<3,3>::Type::BlockType Block ;

	for( int j = 0 ; j < MeshType::NV ; ++j ) {
		for( int k = 0 ; k < MeshType::NV ; ++k ) {
			Block &b = A.block( rowIndices[itp.nodes[k]], colIndices[itp.nodes[j]] ) ;

			// dux_dx, duy_dy, duz_dz
			b(0,0) +=      w * dc_dx(k, 0) * dc_dx(j, 0)
				   +  .5 * w * dc_dx(k, 1) * dc_dx(j, 1)
				   +  .5 * w * dc_dx(k, 2) * dc_dx(j, 2) ;
			b(1,1) +=      w * dc_dx(k, 1) * dc_dx(j, 1)
				   +  .5 * w * dc_dx(k, 0) * dc_dx(j, 0)
				   +  .5 * w * dc_dx(k, 2) * dc_dx(j, 2) ;
			b(2,2) +=      w * dc_dx(k, 2) * dc_dx(j, 2)
				   +  .5 * w * dc_dx(k, 0) * dc_dx(j, 0)
				   +  .5 * w * dc_dx(k, 1) * dc_dx(j, 1) ;

			// 2 * .25 * (dux_dy + duy_dx) * (dvx_dy + dvy_dx)
			b(0,1) += .5 * w * dc_dx(k, 1) * dc_dx(j, 0) ;
			b(1,0) += .5 * w * dc_dx(k, 0) * dc_dx(j, 1) ;

			// 2 * .25 * (dux_dz + duz_dx) * (dvx_dz + dvz_dx)
			b(0,2) += .5 * w * dc_dx(k, 2) * dc_dx(j, 0) ;
			b(2,0) += .5 * w * dc_dx(k, 0) * dc_dx(j, 2) ;

			// 2 * .25 * (duz_dy + duy_dz) * (dvz_dy + dvy_dz)
			b(2,1) += .5 * w * dc_dx(k, 1) * dc_dx(j, 2) ;
			b(1,2) += .5 * w * dc_dx(k, 2) * dc_dx(j, 1) ;

		}
	}
}

} //d6
