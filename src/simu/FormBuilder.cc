#include "FormBuilder.hh"

#include "geo/Grid.hh"
#include "geo/Tensor.hh"

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

void FormBuilder::addTauDu( FormMat<6,3>::Type& A, Scalar w, Itp itp, Dcdx dc_dx, Indices rowIndices, Indices colIndices )
{
	typedef FormMat<6,3>::Type::BlockType Block ;

//#pragma omp parallel for
	for( int k = 0 ; k < MeshType::NV ; ++k ) {
		const Scalar m = w * itp.coeffs[k] ;
		for( int j = 0 ; j < MeshType::NV ; ++j ) {
			Block &b = A.block( rowIndices[itp.nodes[k]], colIndices[itp.nodes[j]] ) ;

			// a * sqrt2_3 * (dux_dx + duy_dy + duz_dz)
			b(0,0) += m * s_sqrt_23 * dc_dx(j, 0) ;
			b(0,1) += m * s_sqrt_23 * dc_dx(j, 1) ;
			b(0,2) += m * s_sqrt_23 * dc_dx(j, 2) ;

			// b * (dux_dx - duy_dy )
			b(1,0) += m * dc_dx(j, 0) ;
			b(1,1) -= m * dc_dx(j, 1) ;

			// c * isqrt_3 * ( -dux_dx - duy_dy + 2*duz_dz)
			b(2,0) -= m * s_isqrt_3 * dc_dx(j, 0) ;
			b(2,1) -= m * s_isqrt_3 * dc_dx(j, 1) ;
			b(2,2) += m * s_isqrt_3 * dc_dx(j, 2) * 2;

			// d * ( dux_dy + duy_dx )
			b(3,0) += m * dc_dx(j, 1) ;
			b(3,1) += m * dc_dx(j, 0) ;

			// e * ( dux_dz + duz_dx )
			b(4,0) += m * dc_dx(j, 2) ;
			b(4,2) += m * dc_dx(j, 0) ;

			// f * ( duz_dy + duy_dz )
			b(5,2) += m * dc_dx(j, 1) ;
			b(5,1) += m * dc_dx(j, 2) ;
		}
	}
}

void FormBuilder::addVDp ( FormMat<3,1>::Type& A, Scalar w, Itp itp, Dcdx dc_dx, Indices rowIndices, Indices colIndices )
{
	typedef FormMat<3,1>::Type::BlockType Block ;

//#pragma omp parallel for
	for( int k = 0 ; k < MeshType::NV ; ++k ) {
		const Scalar m = w * itp.coeffs[k] ;
		for( int j = 0 ; j < MeshType::NV ; ++j ) {
			Block &b = A.block( rowIndices[itp.nodes[k]], colIndices[itp.nodes[j]] ) ;

			// ( vx da_dx + vy da_dy + vz da_dz)
			b += m * dc_dx.row(j) ;
		}
	}
}

void FormBuilder::addTauWu( FormMat<3,3>::Type& A, Scalar w, Itp itp, Dcdx dc_dx, Indices rowIndices, Indices colIndices )
{
	typedef FormMat<3,3>::Type::BlockType Block ;

//#pragma omp parallel for
	for( int k = 0 ; k < MeshType::NV ; ++k ) {
		const Scalar m = w * itp.coeffs[k] ;
		for( int j = 0 ; j < MeshType::NV ; ++j ) {
			Block &b = A.block( rowIndices[itp.nodes[k]], colIndices[itp.nodes[j]] ) ;

			// i * ( dux_dy - duy_dx )
			b(0,0) += m * dc_dx(j, 1) ;
			b(0,1) -= m * dc_dx(j, 0) ;

			// j * ( dux_dz - duz_dx )
			b(1,0) += m * dc_dx(j, 2) ;
			b(1,2) -= m * dc_dx(j, 0) ;

			// k * ( duy_dz - duz_dy )
			b(2,1) += m * dc_dx(j, 2) ;
			b(2,2) -= m * dc_dx(j, 1) ;
		}
	}
}

void FormBuilder::addUTauGphi( FormMat<6,3>::Type& A, Scalar w, Itp itp, const Vec& dphi_dx, Indices rowIndices, Indices colIndices )
{
	typedef FormMat<6,3>::Type::BlockType Block ;

//#pragma omp parallel for
	for( int k = 0 ; k < MeshType::NV ; ++k ) {
		for( int j = 0 ; j < MeshType::NV ; ++j ) {
			Block &b = A.block( rowIndices[itp.nodes[k]], colIndices[itp.nodes[j]] ) ;
			const Scalar m = w * itp.coeffs[k] * itp.coeffs[j] ;

			// a * sqrt2_3 * (dphi_dx ux + dphi_dy uy + dphi_dz uz)
			b(0,0) += m * s_sqrt_23 * dphi_dx(0) ;
			b(0,1) += m * s_sqrt_23 * dphi_dx(1) ;
			b(0,2) += m * s_sqrt_23 * dphi_dx(2) ;

			// b * (dphi_dx ux - dphi_dy uy )
			b(1,0) += m * dphi_dx(0) ;
			b(1,1) -= m * dphi_dx(1) ;

			// c * isqrt_3 * ( -dux_dx - duy_dy + 2*duz_dz)
			b(2,0) -= m * s_isqrt_3 * dphi_dx(0) ;
			b(2,1) -= m * s_isqrt_3 * dphi_dx(1) ;
			b(2,2) += m * s_isqrt_3 * dphi_dx(2) * 2;

			// d * ( dux_dy + duy_dx )
			b(3,0) += m * dphi_dx(1) ;
			b(3,1) += m * dphi_dx(0) ;

			// e * ( dux_dz + duz_dx )
			b(4,0) += m * dphi_dx(2) ;
			b(4,2) += m * dphi_dx(0) ;

			// f * ( duz_dy + duy_dz )
			b(5,2) += m * dphi_dx(1) ;
			b(5,1) += m * dphi_dx(2) ;
		}
	}

}

} //d6
