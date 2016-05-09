#include "DiphasicPrimal.hh"

#include "utils/Log.hh"

#include <bogus/Core/Block.impl.hpp>
#include <bogus/Core/Block.io.hpp>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include <fstream>

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

	const Index m = this->m() ;
	const Index r = this->r() ;
	const Index p = this->p() ;
	const Index s = this->s() ;

	const Index padding = this->padding() ;

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

	M.resize( s, s ) ;

	typedef Eigen::Triplet<Scalar> Tpl ;
	std::vector< Tpl > tpl ;
	tpl.reserve( A.nonZeros() + R.nonZeros() + p + 2*C.nonZeros() + 2*B.nonZeros()  );

	for( Index i = 0 ; i < m ; ++i ) {
		for( SM::InnerIterator it (A, i) ; it ; ++it )
		{
			tpl.push_back( Tpl( it.row(), i, it.value() ) );
		}
		for( SM::InnerIterator it (B, i) ; it ; ++it )
		{
			tpl.push_back( Tpl( m+r + it.row(), i, -it.value() ) );
			tpl.push_back( Tpl( i, m+r + it.row(), -it.value() ) );
		}
	}
	for( Index i = 0 ; i < r ; ++i ) {
		for( SM::InnerIterator it (R, i) ; it ; ++it )
		{
			tpl.push_back( Tpl( m+it.row(), m+i, it.value() ) );
		}
		for( SM::InnerIterator it (C, i) ; it ; ++it )
		{
			tpl.push_back( Tpl( m+r + it.row(), m+i, -it.value() ) );
			tpl.push_back( Tpl( m+i, m+r + it.row(), -it.value() ) );
		}
	}
	for( Index i = 0 ; i < p ; ++i ) {
//		for( SM::InnerIterator it (Q, i) ; it ; ++it )
//		{
//			tpl.push_back( Tpl( m+r+it.row(), m+r+i, -pen * it.value() ) );
//		}
		tpl.push_back( Tpl( m+r + i, m+r + i, - pen ) );
	}

	for( Index i = 0 ; i < padding ; ++i ) {
		tpl.push_back( Tpl( m+r+p+i, m+r+p+i, 1 ) );
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

template <class Archive>
void DiphasicPrimalData::serialize(Archive &ar, const unsigned int )
{
	ar & A ;
	ar & M_lumped_inv ;
	ar & R ;

	ar & B ;
	ar & C ;
	ar & G ;
	ar & H ;

	ar & k ;

	ar & mu ;

}

bool DiphasicPrimalData::load(const char *file)
{
	std::ifstream ifs( file );

	if(!ifs) {
		Log::Error() << "Cannot read " << file << std::endl ;
		return false ;
	}

	boost::archive::binary_iarchive ia(ifs);
	ia >> (*this) ;

	return true ;
}

bool DiphasicPrimalData::dump(const char *file) const
{
	std::ofstream ofs( file );

	if(!ofs) {
		Log::Error() << "Cannot write " << file << std::endl ;
		return false ;
	}

	boost::archive::binary_oarchive oa(ofs);
	oa << (*this) ;
	return true ;
}


} //d6
