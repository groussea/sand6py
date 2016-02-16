#ifndef D6_RIGID_BODY_DATA_HH
#define D6_RIGID_BODY_DATA_HH

#include "ActiveIndices.hh"
#include "PhaseFields.hh"

#include "geo/BoundaryInfo.hh"

#include "utils/block_mat.hh"

namespace d6 {

class RigidBody ;

struct RigidBodyData
{
	typedef RBStresses TensorField ;

	RigidBodyData( RigidBody& rb_, TensorField &s ) ;

	//! Volume fraction taken by the rigid-body at a given point
	Scalar    phi( const Vec &x ) const ;
	//! Gradient of the volume fraction taken by the rigid-body at a given point
	void grad_phi( const Vec &x, Vec &grad ) const ;

	//! Computes nodes that are influenced by the rigid-body
	void compute_active( const Active& phaseNodes ) ;
	//! Assembles projection and jacobian matrices
	void assemble_matrices(const PrimalShape &primalShape, const DualShape &dualShape,
						   const Active &primalNodes, const Active& dualNodes,
						   Index totNodes ) ;

	RigidBody&   rb ;
	TensorField& stresses ;

	Active	    nodes ;

	FormMat<SD,WD>::Type	jacobian ;     //!< int( (u grad phi):tau )
	FormMat<SD,WD>::Type	projection ;   //!< Linear operator giving rb velocities at mesh nodes

	//! Integrate volume fraction at dual nodes
	DynArr intFraction ;
private:
	static const Scalar s_splatRad ;

	void integrate( const PrimalShape& primalShape, const DualShape& dualShape,
					const Active &primalNodes, const Active& dualNodes,
					Index totNodes  ) ;

};

} //d6

#endif
