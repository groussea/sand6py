#ifndef D6_DIPHASIC_STEP_DATA_HH
#define D6_DIPHASIC_STEP_DATA_HH

#include "mono/PhaseFields.hh"
#include "mono/PhaseStepData.hh"

#include "simu/ActiveIndices.hh"
#include "simu/FormBuilder.hh"
#include "simu/StressMassMatrix.hh"

namespace d6 {

struct Config ;
class DynParticles ;
struct Phase ;

struct DiphasicStepData {

	//! Active nodes: nodes under the influence of at least one particle
	Active primalNodes ;
	Active dualNodes ;

	struct Forms {
		// Linear form vectors

		DynVec linearMomentum ;  //!< integral of fraction times velocity

		DynVec externalForces ; //!< integral of external forces
		DynArr volumes  ;	    //!< Volume aossciated to each node = int(1)

		DynArr fraction  ;		//!< integrated volume fraction occupied by the phase

		// Bilinear Form marices

		//! Lumped mass matrix and its inverse
		typename FormMat<WD,WD>::SymType M_lumped ;
		typename FormMat<WD,WD>::SymType M_lumped_inv ;

		typename FormMat<WD,WD>::Type A ; //!< Mass + Visco ; Could be Symmetric when FormBuilder has sym index

		typename FormMat< 1,WD>::Type B ; //!< p div v
		typename FormMat< 1,WD>::Type C ; //!< grad(p) w alpha phi

		typename FormMat<SD,WD>::Type G ; //!< \phi Tau:D(u)
		typename FormMat<SD,WD>::Type H ; //!< \phi Tau:D(w)


		typedef AbstractStressMassMatrix< DualShape > StressMassMatrix ;
		StressMassMatrix S ;

	} forms ;

	//! Projectors enforcing boundary conditions
	typedef PhaseStepData::Projectors Projectors ;
	Projectors activeProj ;
	Projectors fullGridProj ;

	Index nPrimalNodes() const
	{
		return primalNodes.count() ;
	}

	Index nDualNodes() const
	{
		return dualNodes.count() ;
	}

	void compute(
			const DynParticles& particles, const Config &config, const Scalar dt,
			Phase &phase  ) ;


private:

	void assembleMatrices(const Particles& particles,
		const Config &config, const Scalar dt, const DualShape &dShape, const Phase& phase,
		const PrimalScalarField &intPhi, const PrimalVectorField &intPhiVel ) ;

	static void computeProjectors(const Config &config, const PrimalShape &pShape,
						   Projectors& mats ) ;
};

} //d6


#endif
