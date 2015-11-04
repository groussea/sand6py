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
	void read(std::vector< bool > &activeCells,
			   ScalarField &phi, VectorField &phiVel,
			   ScalarField &phiInertia, TensorField &phiOrient
			   ) const ;

	const Particles &geo() const { return m_geo ; }
	Particles &geo() { return m_geo ; }

	size_t count() const { return m_geo.count() ; }

	void clamp_particle( size_t i, const MeshType &mesh ) ;

private:

	void splitMerge( const MeshType & mesh ) ;

	void resize( size_t n ) ;

	Particles m_geo ;
	Scalar m_meanVolume ;

	Eigen::Matrix< Scalar, 9, Eigen::Dynamic > m_affine ;
	Eigen::Matrix< Scalar, 1, Eigen::Dynamic > m_inertia ;

};

} //d6

#endif
