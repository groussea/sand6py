#ifndef D6_DYN_PARTICLES_HH
#define D6_DYN_PARTICLES_HH

#include "geo/Particles.hh"

class PhaseDescr ;

namespace d6 {

class DynParticles {

public:
	DynParticles() ;

	void generate( const Config &c, const MeshType& mesh ) ;

	const Particles &geo() const { return m_geo ; }
	size_t count() const { return m_geo.count() ; }

private:

	void resize( size_t n ) ;

	Particles m_geo ;


	Eigen::Matrix< Scalar, 9, Eigen::Dynamic > m_affine ;

};

} //d6

#endif
