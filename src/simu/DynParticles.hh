#ifndef D6_DYN_PARTICLES_HH
#define D6_DYN_PARTICLES_HH

#include "simu/PhaseFields.hh"

#include "geo/Particles.hh"

namespace d6 {

struct Config ;
struct Phase ;
class Scenario ;

class DynParticles {

public:
	DynParticles() ;

	void generate( const Config &c, const MeshType& mesh, const Scenario &scenario ) ;

	void update( const Config&c, const Scalar dt, const Phase& phase ) ;

	void integratePrimal( std::vector< bool > &activeCells,
			   PrimalScalarField &phi,    PrimalVectorField &phiVel ) const ;

	void integrateDual( std::vector< bool > &activeCells,
			   DualScalarField &phi,      DualScalarField &phiInertia,
			   DualTensorField &phiOrient,DualScalarField &phiCohesion
			   ) const ;

	const Particles &geo() const { return m_geo ; }
	Particles &geo() { return m_geo ; }

	size_t count() const { return m_geo.count() ; }

	void remove( size_t particleId ) ;

	const Particles::EventLog& events() const {
		return m_events ;
	}
	Particles::EventLog& events() {
		return m_events ;
	}

	const Particles::Data< 1 >::Type&  inertia() const { return m_inertia ; }
	const Particles::Data< 1 >::Type& cohesion() const { return m_cohesion ; }

private:

	void clamp_particle( size_t i, const MeshType &mesh ) ;

	void splitMerge( const MeshType & mesh ) ;

	void resize( size_t n ) ;

	Particles m_geo ;
	Particles::EventLog m_events ;
	Scalar m_meanVolume ;

	Particles::Data< WD*WD >::Type m_affine ;
	Particles::Data< 1 >::Type m_inertia ;
	Particles::Data< 1 >::Type m_cohesion ;

};

} //d6

#endif
