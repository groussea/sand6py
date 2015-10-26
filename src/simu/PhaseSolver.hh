#ifndef D6_PHASE_SOLVE_HH
#define D6_PHASE_SOLVE_HH

#include "geo/geo.fwd.hh"

namespace d6 {

class DynParticles ;
class Config ;
struct Phase ;

class PhaseSolver {

public:
	explicit PhaseSolver( const DynParticles& particles ) ;

	void step(const Config &config, const MeshType& mesh, Phase& phase ) ;

private:

	const DynParticles& m_particles ;

};


} //d6


#endif
