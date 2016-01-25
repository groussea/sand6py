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

	//! Volume fraction taken by the rigid-body at a given point
	Scalar    phi( const Vec &x ) const ;
	//! Gradient of the volume fraction taken by the rigid-body at a given point
	void grad_phi( const Vec &x, Vec &grad ) const ;

	//! Computes nodes that are influenced by the rigid-body
	void compute_active( const Active& phaseNodes, BoundaryConditions &bc ) ;
	//! Assembles projection and jacobian matrices
	void assemble_matrices( const Active& phaseNodes, Index totNodes ) ;

	RigidBody&   rb ;
	TensorField& stresses ;

	Active	    nodes ;
	typename MeshType::Cells occupiedCells ;

	FormMat<SD,WD>::Type	jacobian ;     //!< int( (u grad phi):tau )
	FormMat<SD,WD>::Type	projection ;   //!< Linear operator giving rb velocities at mesh nodes

	DynVec fraction ;  //!< Interpolated volume fraction at occupied nodes

private:
	static const Scalar s_splatRad ;

	void integrate(const Active& phaseNodes, Index totNodes  ) ;

};

} //d6

#endif
