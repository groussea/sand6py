#ifndef D6_PHASE_SOLVE_HH
#define D6_PHASE_SOLVE_HH

#include "ActiveIndices.hh"

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

	~PhaseSolver() ;

	void step(const Config &config, const Scalar dt,
			  Phase& phase, Stats &stats,
			  std::vector<RigidBody> &rigidBodies,
			  std::vector<TensorField> &rbStresses
			  ) const ;

private:



	void solve(const Config& config, const Scalar dt,
			   const PhaseStepData &stepData,
			   Phase& phase, std::vector< RigidBodyData > &rbData, Stats& stats ) const ;

	void addRigidBodyContrib(const Config &c, const Scalar dt, const PhaseStepData &stepData,
							 const DynVec &u, const RigidBodyData &rb,
							 PrimalData& primalData, DynArr &totFraction ) const ;
	void addCohesionContrib (const Config&c, const PhaseStepData &stepData,
							  PrimalData& primalData, DynVec &u ) const ;

	void solveComplementarity(const Config&c, const Scalar dt,
							  const PhaseStepData& stepData ,
							  std::vector< RigidBodyData >& rbData,
							  DynVec &u, Phase &phase, Stats &stats ) const ;

	void enforceMaxFrac(const Config &c, const PhaseStepData &stepData,
									   const std::vector<RigidBodyData> &rbData,
									   DynVec &depl ) const ;


	const DynParticles& m_particles ;

};


} //d6


#endif
