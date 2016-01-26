#include "FormBuilder.hh"

#include "geo/MeshImpl.hh"

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

void FormBuilder::reset(Index rows)
{
	m_data.clear();
	m_data.resize( rows );
	m_compressed.outer.clear() ;
}

void FormBuilder::addRows( Index rows )
{
	m_data.resize( m_data.size() + rows );

	const Index old = m_compressed.outer.size() ;
	if( old > 0  )
	{
		m_compressed.outer.resize( old + rows ) ;
		for( Index i = 0 ; i < rows ; ++i  ) {
			m_compressed.outer[ old+i ] = m_compressed.outer[ old-1 ] ;
		}
	}
}

void FormBuilder::addDuDv( FormMat<WD,WD>::Type& A, Scalar w, Itp itp, Dcdx dc_dx, Indices rowIndices, Indices colIndices )
{
	for( int k = 0 ; k < MeshType::NV ; ++k ) {
		addDuDv( A, w, rowIndices[itp.nodes[k]], dc_dx.row(k), itp, dc_dx, colIndices );
	}
}

void FormBuilder::addTauDu( FormMat<SD,WD>::Type& A, Scalar w, Itp itp, Dcdx dc_dx, Indices rowIndices, Indices colIndices )
{
	for( int k = 0 ; k < MeshType::NV ; ++k ) {
		addTauDu( A, w * itp.coeffs[k], rowIndices[itp.nodes[k]], itp, dc_dx, colIndices );
	}
}

void FormBuilder::addVDp( FormMat<WD,1>::Type& A, Scalar w, Itp itp, Dcdx dc_dx, Indices rowIndices, Indices colIndices )
{
	for( int k = 0 ; k < MeshType::NV ; ++k ) {
		addVDp( A, w * itp.coeffs[k], rowIndices[itp.nodes[k]], itp, dc_dx, colIndices );
	}
}

void FormBuilder::addTauWu( FormMat<RD,WD>::Type& A, Scalar w, Itp itp, Dcdx dc_dx, Indices rowIndices, Indices colIndices )
{
	for( int k = 0 ; k < MeshType::NV ; ++k ) {
		addTauWu( A, w * itp.coeffs[k], rowIndices[itp.nodes[k]], itp, dc_dx, colIndices );
	}
}

} //d6
