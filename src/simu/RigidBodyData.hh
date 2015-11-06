#ifndef D6_RIGID_BODY_DATA_HH
#define D6_RIGID_BODY_DATA_HH

#include "utils/block_mat.hh"
#include "ActiveIndices.hh"

#include "geo/BoundaryInfo.hh"

namespace d6 {

class RigidBody ;

struct RigidBodyData
{
	RigidBodyData( RigidBody& rb_, TensorField &s ) ;

	Scalar phi( const Vec &x ) const ;

	void compute_active( const Active& phaseNodes, BoundaryConditions &bc ) ;

	RigidBody&   rb ;
	TensorField& stresses ;

	Active	    nodes ;

	FormMat<6,3>::Type	m_jacobian ;
	FormMat<6,3>::Type	m_projection ;

};

} //d6

#endif
