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


} //d6
