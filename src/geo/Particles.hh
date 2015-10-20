#ifndef D6_PARTICLES_HH
#define D6_PARTICLES_HH

#include "geo.fwd.hh"
#include "utils/alg.hh"

namespace d6 {

class Config ;

class Particles {

public:

	static const size_t s_MAX ;

	Particles() ;

	void generate( const Config &c, const MeshType& mesh ) ;

	size_t count() const { return m_count ; }

private:

	std::size_t m_count ;

	Eigen::VectorXd m_volumes ;

	Eigen::Matrix< Scalar, 3, Eigen::Dynamic > m_centers ;
	Eigen::Matrix< Scalar, 3, Eigen::Dynamic > m_velocities ;

	Eigen::Matrix< Scalar, 6, Eigen::Dynamic > m_frames ;

	Eigen::Matrix< Scalar, 6, Eigen::Dynamic > m_orient ; // Aniso


	void resize( size_t n ) ;
};

} //d6

#endif

