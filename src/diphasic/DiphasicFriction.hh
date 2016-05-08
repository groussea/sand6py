#ifndef D6_DIPHASIC_FRICTION_HH
#define D6_DIPHASIC_FRICTION_HH

#include "utils/alg.hh"
#include "utils/block_mat.hh"

#include <Eigen/Sparse>

namespace d6 {

struct DiphasicPrimalData {
	typedef FormMat< SD, WD >::Type HType ;

	HType H ;
	HType G ;

	DynVec k ;

	DynVec mu ;

};

struct DiphasicFrictionSolver {

	typedef Eigen::SparseMatrix< Scalar > ESM ;
	typedef Eigen::SimplicialLDLT< ESM > ELDLT ;

	typedef DiphasicPrimalData::HType HType ;

	// x = u w p
	// M (u w) = l + (G H)' lambda
	// gamma = Gu + Hw + k
	// (lambda, gamma) \in DPmu

	Scalar solve(const ESM &M,
				  const DiphasicPrimalData& data,
				  DynVec &x, DynVec &lambda
				) ;


};


} //d6

#endif
