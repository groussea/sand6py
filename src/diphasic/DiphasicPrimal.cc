#include "DiphasicPrimal.hh"

#include <bogus/Core/Block.io.hpp>

#include <bogus/Core/Block.impl.hpp>

namespace d6 {

void DiphasicPrimalData::makePenalizedEigenStokesMatrix(
		Eigen::SparseMatrix<Scalar> & M,
		const Scalar pen
		) const
{
	// Assemble Eigen mat for
	// (  A    0     -B^T  )
	// (  0    R     -C^T  )
	// ( -B   -C  -pen Id  )

	const Index m  = A.rows() ;
	const Index ma = R.rows() ;
	const Index n  = B.rows() ;
	const Index r = m + n + ma;

	typedef Eigen::SparseMatrix< Scalar > SM ;
	SM A, B, C, R, Q ;
	bogus::convert( this->A, A ) ;
	bogus::convert( this->R, R ) ;

	bogus::convert( this->B, B ) ;
	bogus::convert( this->C, C ) ;


	A.prune(1.) ;
	B.prune(1.) ;
	C.prune(1.) ;
	R.prune(1.) ;

	M.resize( r, r ) ;

	typedef Eigen::Triplet<Scalar> Tpl ;
	std::vector< Tpl > tpl ;
	tpl.reserve( A.nonZeros() + R.nonZeros() + n + 2*C.nonZeros() + 2*B.nonZeros()  );

	for( Index i = 0 ; i < m ; ++i ) {
		for( SM::InnerIterator it (A, i) ; it ; ++it )
		{
			tpl.push_back( Tpl( it.row(), i, it.value() ) );
		}
		for( SM::InnerIterator it (B, i) ; it ; ++it )
		{
			tpl.push_back( Tpl( m+ma + it.row(), i, -it.value() ) );
			tpl.push_back( Tpl( i, m+ma + it.row(), -it.value() ) );
		}
	}
	for( Index i = 0 ; i < ma ; ++i ) {
		for( SM::InnerIterator it (R, i) ; it ; ++it )
		{
			tpl.push_back( Tpl( m+it.row(), m+i, it.value() ) );
		}
		for( SM::InnerIterator it (C, i) ; it ; ++it )
		{
			tpl.push_back( Tpl( m+ma + it.row(), m+i, -it.value() ) );
			tpl.push_back( Tpl( m+i, m+ma + it.row(), -it.value() ) );
		}
	}
	for( Index i = 0 ; i < n ; ++i ) {
//		for( SM::InnerIterator it (Q, i) ; it ; ++it )
//		{
//			tpl.push_back( Tpl( m+ma+it.row(), m+ma+i, -pen * it.value() ) );
//		}
		tpl.push_back( Tpl( m+ma + i, m+ma + i, - pen ) );
	}

	M.setFromTriplets( tpl.begin(), tpl.end() ) ;

}

void  DiphasicPrimalData::factorize( ESM &mat, MInvType &sbm )
{

	std::vector< unsigned > rpb {(unsigned)  mat.rows()} ;
	std::vector< unsigned > cpb {(unsigned)  mat.cols()} ;

	sbm.setRows( rpb );
	sbm.setCols( cpb );

	sbm.insertBack(0,0).compute( mat ) ;
	sbm.finalize();

}


} //d6
