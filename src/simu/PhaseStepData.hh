#ifndef D6_PHASE_STEP_DATA_HH
#define D6_PHASE_STEP_DATA_HH

#include "ActiveIndices.hh"
#include "FormBuilder.hh"
#include "RigidBodyData.hh"

#include "PhaseFields.hh"
#include "StressMassMatrix.hh"

#include "geo/BoundaryInfo.hh"

namespace d6 {

struct Config ;
class RigidBody ;
class DynParticles ;
struct Phase ;

struct PhaseStepData {

	//! Active nodes: nodes under the influence of at least one particle
	Active primalNodes ;
	Active dualNodes ;

	struct Forms {
		// Linear form vectors

		DynVec phiu ;           //!< integral of fraction times velocity
		DynVec externalForces ; //!< integral of external forces
		DynArr volumes  ;	    //!< Volume aossciated to each node = int(1)

		DynArr fraction  ;		//!< integrated volume fraction occupied by the phase

		// Bilinear Form marices

		//! Lumped mass matrix, its inverse and a factorization
		typename FormMat<WD,WD>::SymType M_lumped ;
		typename FormMat<WD,WD>::SymType M_lumped_inv ;
		typename FormMat<WD,WD>::SymType M_lumped_inv_sqrt ;

		typename FormMat<WD,WD>::Type A ; //!< Mass + Visco ; Could be Symmetric when FormBuilder has sym index

		typename FormMat<SD,WD>::Type B ; //!< \phi Tau:D(u)
		typename FormMat<RD,WD>::Type J ; //!< \phi Tau:W(u)

		typename FormMat<1 ,WD>::Type C ; //!< \phi v.grad(p)

		StressMassMatrix S ;

	} forms ;

	//! Projectors enforcing boundary conditions
	struct Projectors {
		typename FormMat<WD,WD>::SymType vel ;
		typename FormMat<SD,SD>::SymType stress ;
	} proj ;


	DynArr cohesion ;   //!< interpolated cohesivity
	DynArr inertia  ;   //!< interpolated inertial number

	typename FormMat<SD,SD>::SymType Aniso ; //!< Anisotropy linear operator (N^{-1})

	PhaseStepData()
		: m_totRbNodes( 0 )
	{}

	Index nPrimalNodes() const
	{
		return primalNodes.count() ;
	}

	Index nDualNodes() const
	{
		return dualNodes.count() ;
	}

	Index nSuppNodes() const
	{
		return m_totRbNodes ;
	}

	void compute(
			const DynParticles& particles, const Config &config, const Scalar dt,
			Phase &phase,
			std::vector< RigidBody   >& rigidBodies,
			std::vector< RBStresses > &rbStresses,
			std::vector< RigidBodyData > &rbData  ) ;

private:
	void computeActiveNodes(const std::vector< bool >& activeCells,
							const PrimalShape &pShape , const DualShape &dShape) ;
	void computeActiveBodies( std::vector<RigidBody> &rigidBodies,
							  std::vector<RBStresses> &rbStresses,
							  std::vector< RigidBodyData > &rbData ) ;

	void computeProjectors(const Config &config, const PrimalShape &pShape, const DualShape &dShape,
						   const std::vector<RigidBodyData> & rbData, Projectors& mats ) const ;

	void computeAnisotropy(const DynVec& orientation,  const Config &config,
						   typename FormMat<SD,SD>::SymType &Aniso ) const ;

	void computePhiAndGradPhi(const PrimalScalarField& intPhi, PrimalScalarField&phi, PrimalVectorField &grad_phi ) const ;

	void assembleMatrices(const Particles& particles, const Config& c, const Scalar dt,
						   const DualShape &dShape, const PrimalScalarField &phiInt,
						   std::vector< RigidBodyData > &rbData ) ;


	PhaseStepData( const PhaseStepData& ) = delete ;
	PhaseStepData& operator=( const PhaseStepData& ) = delete ;

	Index m_totRbNodes ;
};


} //d6


#endif
