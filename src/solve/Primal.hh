#ifndef D6_PRIMAL_HH
#define D6_PRIMAL_HH

#include "utils/block_mat.hh"

namespace  d6 {

struct PrimalData {
	typedef typename FormMat<6,3>::Type HType ;

	HType H   ;
	DynVec w  ;

	DynVec mu ;

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
