#ifndef D6_DIPHASIC_SOLVER_HH
#define D6_DIPHASIC_SOLVER_HH

#include "utils/scalar.hh"

namespace d6 {

class DynParticles ;

struct Phase ;
struct Config ;
struct DiphasicStepData ; ;

class DiphasicSolver {


public:
	explicit DiphasicSolver(
			const DynParticles& particles
			) ;


	//! Solve for end-of-steps velocities, reading initial quantities from the particles
	void step(const Config &config, const Scalar dt,
			  Phase& phase  ) const ;


private:
	const DynParticles& m_particles ;

	void solve(
			const Config& config, const Scalar dt, const DiphasicStepData& stepData ,
			Phase& phase ) const ;
};


} //d6

#endif
