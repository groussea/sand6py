#ifndef D6_DYN_PARTICLES_HH
#define D6_DYN_PARTICLES_HH

#include "geo/Particles.hh"

namespace d6 {

class Config ;
struct Phase ;

class DynParticles {

public:
	DynParticles() ;

	void generate( const Config &c, const MeshType& mesh ) ;

	void update( const Config&c, const Phase& phase ) ;

	const Particles &geo() const { return m_geo ; }
	size_t count() const { return m_geo.count() ; }

	void clamp_particle( size_t i, const MeshType &mesh ) ;

private:

	void resize( size_t n ) ;

	Particles m_geo ;


	Eigen::Matrix< Scalar, 9, Eigen::Dynamic > m_affine ;
	Eigen::Matrix< Scalar, 1, Eigen::Dynamic > m_inertia ;

};

} //d6

#endif
