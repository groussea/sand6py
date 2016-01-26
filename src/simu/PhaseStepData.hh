#ifndef D6_PHASE_STEP_DATA_HH
#define D6_PHASE_STEP_DATA_HH

#include "ActiveIndices.hh"
#include "FormBuilder.hh"
#include "RigidBodyData.hh"

#include "geo/BoundaryInfo.hh"

namespace d6 {

struct Config ;
class RigidBody ;
class DynParticles ;
struct Phase ;

struct PhaseStepData {

	//! Active nodes: nodes under the influence of at least one particle
	Active nodes ;

	struct Forms {
		// Linear form vectors

		DynVec phiu ;           //!< integral of fraction times velocity
		DynVec externalForces ; //!< integral of external forces
		DynVec volumes  ;	    //!< Volume aossciated to each node = int(1)

		// Bilinear Form marices

		//! Lumped mass matrix, its inverse and a factorization
		typename FormMat<WD,WD>::SymType M_lumped ;
		typename FormMat<WD,WD>::SymType M_lumped_inv ;
		typename FormMat<WD,WD>::SymType M_lumped_inv_sqrt ;

		typename FormMat<WD,WD>::Type A ; //!< Mass + Visco ; Could be Symmetric when FormBuilder has sym index

		typename FormMat<SD,WD>::Type B ; //!< \phi Tau:D(u)
		typename FormMat<RD,WD>::Type J ; //!< \phi Tau:W(u)

		typename FormMat<WD,1>::Type C ; //!< \phi v.grad(p)
	} forms ;

	//! Projectors enforcing boundary conditions
	struct Projectors {
		typename FormMat<WD,WD>::SymType vel ;
		typename FormMat<SD,SD>::SymType stress ;
	} proj ;


	DynVec fraction  ;	//!< interpolated volume fraction occupied by the phase
	DynVec cohesion ;   //!< interpolated cohesivity
	DynVec inertia  ;   //!< interpolated inertial number

	typename FormMat<SD,SD>::SymType Aniso ; //!< Anisotropy linear operator (N^{-1})

	PhaseStepData()
		: m_totRbNodes( 0 )
	{}

	Index nNodes() const
	{
		return nodes.count() ;
	}

	Index nSuppNodes() const
	{
		return m_totRbNodes ;
	}

	void compute(
			const DynParticles& particles, const Config &config, const Scalar dt,
			Phase &phase,
			std::vector< RigidBody   >& rigidBodies,
			std::vector<TensorField > &rbStresses,
			std::vector< RigidBodyData > &rbData  ) ;

private:
	void computeActiveNodes(const std::vector< bool >& activeCells,
							const VectorField &grad_phi ) ;
	void computeActiveBodies( std::vector<RigidBody> &rigidBodies,
							  std::vector<TensorField> &rbStresses,
							  std::vector< RigidBodyData > &rbData ) ;

	void computeProjectors( const bool weakStressBC, const std::vector<RigidBodyData> & rbData,
							Projectors& mats ) const ;

	void computeAnisotropy(const DynVec& orientation,  const Config &config,
						   typename FormMat<SD,SD>::SymType &Aniso ) const ;

	void computeGradPhi( const ScalarField& fraction, const ScalarField& volumes, VectorField &grad_phi ) const ;

	void assembleMatrices( const Particles& particles, const Config& c, const Scalar dt,
						   const MeshType& mesh, const ScalarField &phiInt,
						   std::vector< RigidBodyData > &rbData ) ;


	PhaseStepData( const PhaseStepData& ) = delete ;
	PhaseStepData& operator=( const PhaseStepData& ) = delete ;

	BoundaryConditions m_surfaceNodes ;
	Index m_totRbNodes ;
};


} //d6


#endif
