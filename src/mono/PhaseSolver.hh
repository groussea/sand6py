#ifndef D6_PHASE_SOLVE_HH
#define D6_PHASE_SOLVE_HH

#include "PhaseFields.hh"

#include "simu/ActiveIndices.hh"

#include "geo/BoundaryInfo.hh"
#include "geo/MeshBase.hh"

#include "utils/scalar.hh"

#include <vector>

namespace d6 {

class DynParticles ;
struct Phase ;
struct PhaseStepData ;

class RigidBody ;
struct RigidBodyData ;
struct PrimalData ;

struct Config ;
class Stats ;

class PhaseSolver {

public:
	explicit PhaseSolver(
			const DynParticles& particles
			) ;

	//! Solve for end-of-steps velocities, reading initial quantities from the particles
	void step(const Config &config, const Scalar dt,
			  Phase& phase, Stats &stats,
			  std::vector<RigidBody> &rigidBodies,
			  std::vector<RBStresses> &rbStresses
			  ) const ;

	static void getCohesiveStress(
			const Config &config, const DynArr& cohesion,  const DynArr& fraction,
			DynVec& cohe_stress ) ;

private:
	//! Two-steps solve of the momentum balance w/ frictional rheology
	void solve(const Config& config, const Scalar dt,
			   const PhaseStepData &stepData,
			   Phase& phase, std::vector< RigidBodyData > &rbData, Stats& stats ) const ;

	//! Assemble and solve the friction problem
	void solveComplementarity(const Config&c, const Scalar dt,
							  const PhaseStepData& stepData ,
							  std::vector< RigidBodyData >& rbData,
							  DynVec &u, Phase &phase, Stats &stats ) const ;

	//! Add contbutions from a rigid body to the friction problem
	void addRigidBodyContrib(const Config &c, const Scalar dt, const PhaseStepData &stepData,
							 const DynVec &u, const RigidBodyData &rb,
							 PrimalData& primalData, DynArr &rbIntFraction ) const ;
	//! Add contribution from chesive forces to the friction problem
	void addCohesionContrib (const Config&c, const PhaseStepData &stepData,
							  PrimalData& primalData, DynVec &u ) const ;

	//! Compute displacement for volume correction
	void enforceMaxFrac(const Config &c, const PhaseStepData &stepData,
									   const std::vector<RigidBodyData> &rbData,
									   DynVec &depl ) const ;


	const DynParticles& m_particles ;

};


} //d6


#endif
