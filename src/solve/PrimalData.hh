#ifndef D6_PRIMAL_DATA_HH
#define D6_PRIMAL_DATA_HH

#include "utils/block_mat.hh"

namespace  d6 {

/*!
 *  Problem defined as
 *  v = H' r
 *  v_rb[j] = inv_im[j] jac[j]' r   \forall j
 *  u = H  v + w + \sum_j{ jac[j] * v_rb[j] }
 *
 *  (u_i, r_i) \in C(mu[i])
 */
struct PrimalData {
	typedef typename FormMat<6,3>::Type HType ;
	typedef typename FormMat<6,6>::Type JacobianType ;
	typedef typename FormMat<6,6>::SymType InvInertiaType ;

	HType H   ;
	DynVec w  ;

	DynVec mu ;

	std::vector< JacobianType > jacobians ;
	std::vector< InvInertiaType > inv_inertia_matrices ;

	Index n() const { return H.rowsOfBlocks() ; }

	template <typename Archive>
	void serialize(Archive &ar, const unsigned int ) ;

	bool load( const char * file ) ;
	bool dump( const char * file ) const ;
};


} //d6

#endif
