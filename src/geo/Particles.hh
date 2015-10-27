#ifndef D6_PARTICLES_HH
#define D6_PARTICLES_HH

#include "Expr.hh"

namespace d6 {

class Particles {

public:

	template< int D >
	struct Data {
		typedef Eigen::Matrix< Scalar, D, Eigen::Dynamic > Type ;
	};


	static const size_t s_MAX ;

	Particles() ;

	void generate(const ScalarExpr &expr, const unsigned nSamples, const MeshType& mesh ) ;

	size_t count() const { return m_count ; }

	template < typename Archive >
	void serialize( Archive &ar, unsigned int ) ;

	const typename Data< 1 >::Type&    volumes() const { return    m_volumes ; }
	const typename Data< 3 >::Type&    centers() const { return    m_centers ; }
	const typename Data< 3 >::Type& velocities() const { return m_velocities ; }
	const typename Data< 6 >::Type&     frames() const { return     m_frames ; }
	const typename Data< 6 >::Type&     orient() const { return     m_orient ; }

private:

	std::size_t m_count ;

	typename Data< 1 >::Type m_volumes ;

	typename Data< 3 >::Type m_centers ;
	typename Data< 3 >::Type m_velocities ;

	typename Data< 6 >::Type m_frames ;
	typename Data< 6 >::Type m_orient ; // Aniso


	void resize( size_t n ) ;

	friend class DynParticles ;

};

} //d6

#endif

