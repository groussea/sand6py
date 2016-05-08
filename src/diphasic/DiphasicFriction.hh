#ifndef D6_DIPHASIC_FRICTION_HH
#define D6_DIPHASIC_FRICTION_HH

#include "DiphasicPrimal.hh"

namespace d6 {


struct DiphasicFrictionSolver {

	typedef Eigen::SparseMatrix< Scalar > ESM ;

	explicit DiphasicFrictionSolver( const DiphasicPrimalData& data)
		:m_data(data)
	{}

	// x = u w p
	// M (u w) = l + (G H)' lambda
	// gamma = Gu + Hw + k
	// (lambda, gamma) \in DPmu

	Scalar solve( const ESM &M, const DiphasicPrimalData::MInvType& M_inv,
				  DynVec &x, DynVec &lambda
				) ;


private:
	const DiphasicPrimalData& m_data ;

};


} //d6

#endif
