#ifndef D6_DIPHASIC_PRIMAL_HH
#define D6_DIPHASIC_PRIMAL_HH

#include "utils/block_mat.hh"

#include <Eigen/Sparse>

namespace d6 {

/*
 *
*/

struct DiphasicPrimalData {

	typedef FormMat< WD, WD >::Type AType ;
	AType A ;
	AType R ;
	AType M_lumped_inv ;

	typedef FormMat<  1, WD >::Type CType ;

	CType B ;
	CType C ;

	typedef FormMat< SD, WD >::Type HType ;

	HType G ;
	HType H ;

	DynVec k ;

	DynVec mu ;

	Index n() const { return H.rowsOfBlocks() ; }
	Index m() const { return A.rows() ; }
	Index r() const { return R.rows() ; }
	Index p() const { return B.rows() ; }
	Index s() const { return m()+r()+p()+padding() ; }

	Index padding() const {
		return (WD - (p()%WD))%WD ;
	}

	template <typename Archive>
	void serialize(Archive &ar, const unsigned int ) ;

	bool load( const char * file ) ;
	bool dump( const char * file ) const ;

	// Not saved


	typedef Eigen::SparseMatrix< Scalar > ESM ;
	typedef bogus::SparseBlockMatrix< bogus::SparseLDLT< Scalar > > MInvType ;

	void makePenalizedEigenStokesMatrix( ESM &M, Scalar pen ) const ;
	static void factorize( ESM &M, MInvType &fact ) ;
};


} //d6

#endif
