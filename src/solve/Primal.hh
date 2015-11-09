#ifndef D6_PRIMAL_HH
#define D6_PRIMAL_HH

#include "utils/block_mat.hh"

namespace  d6 {

struct PrimalData {
	typedef typename FormMat<6,3>::Type HType ;
	typedef typename FormMat<6,6>::Type JacobianType ;

	HType H   ;
	DynVec w  ;

	DynVec mu ;

	std::vector< JacobianType > jacobians ;
	DynMat6 inv_inertia_matrices ;

	Index n() const { return H.rowsOfBlocks() ; }

	bool load( const char * file ) ;
	bool dump( const char * file ) const ;
};

class Primal {

public:
	Primal( const PrimalData &data ) ;

	Scalar solve( DynVec& lambda, DynVec &gamma ) const ;

private:
	const PrimalData& m_data ;

};

} //d6

#endif
