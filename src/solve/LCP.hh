#ifndef D6_LCP_HH
#define D6_LCP_HH

#include "utils/block_mat.hh"

namespace  d6 {

struct LCPData {
	typedef typename FormMat<1,WD>::Type HType ;

	HType H   ;
	DynVec w  ;

	Index n() const { return H.colsOfBlocks() ; }
};

class LCP {

public:
	LCP( const LCPData &data ) ;

	Scalar solve( DynVec& x ) const ;

private:
	const LCPData& m_data ;

};

} //d6

#endif

