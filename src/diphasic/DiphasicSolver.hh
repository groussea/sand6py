#ifndef D6_DIPHASIC_SOLVER_HH
#define D6_DIPHASIC_SOLVER_HH

#include "utils/alg.hh"

namespace d6 {

class DynParticles ;

struct Phase ;
struct FluidPhase ;
struct Config ;
struct DiphasicStepData ;
struct DiphasicPrimalData ;

class DiphasicSolver {


public:
	explicit DiphasicSolver(
			const DynParticles& particles
			) ;


	//! Solve for end-of-steps velocities, reading initial quantities from the particles
	void step(const Config &config, const Scalar dt,
			  Phase& phase, FluidPhase &fluid) const ;


private:
	const DynParticles& m_particles ;

	void solve(const Config& config, const Scalar dt, const DiphasicStepData& stepData ,
			Phase& phase , FluidPhase &fluid) const ;

	void solveStokes( const DiphasicStepData& stepData, const DynVec &l, DynVec &u, DynVec &p ) const ;

	void addCohesionContrib (const Config&c, const DiphasicStepData &stepData,
							DiphasicPrimalData& pbData, DynVec &l ) const ;

};


} //d6

#endif
