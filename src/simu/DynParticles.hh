#ifndef D6_DYN_PARTICLES_HH
#define D6_DYN_PARTICLES_HH

#include "geo/Particles.hh"

namespace d6 {

struct Config ;
struct Phase ;

class DynParticles {

public:
	DynParticles() ;

	void generate( const Config &c, const MeshType& mesh ) ;

	void update( const Config&c, const Phase& phase ) ;
	void read(std::vector< bool > &activeCells,
			   ScalarField &phi, VectorField &phiVel,
			   ScalarField &phiInertia, TensorField &phiOrient,
			   ScalarField &phiCohesion
			   ) const ;

	const Particles &geo() const { return m_geo ; }
	Particles &geo() { return m_geo ; }

	size_t count() const { return m_geo.count() ; }

	void clamp_particle( size_t i, const MeshType &mesh ) ;

	void splitMerge( const MeshType & mesh ) ;

	const Particles::EventLog& events() const {
		return m_events ;
	}
	Particles::EventLog& events() {
		return m_events ;
	}

	const Particles::Data< 1 >::Type&  inertia() const { return m_inertia ; }
	const Particles::Data< 1 >::Type& cohesion() const { return m_cohesion ; }

private:

	void resize( size_t n ) ;

	Particles m_geo ;
	Particles::EventLog m_events ;
	Scalar m_meanVolume ;

	Particles::Data< 9 >::Type m_affine ;
	Particles::Data< 1 >::Type m_inertia ;
	Particles::Data< 1 >::Type m_cohesion ;

};

} //d6

#endif
