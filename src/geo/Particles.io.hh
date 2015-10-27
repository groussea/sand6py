#ifndef D6_PARTICLES_IO_HH
#define D6_PARTICLES_IO_HH

#include "Particles.hh"

#include <boost/serialization/array.hpp>

namespace d6 {

template < typename Archive >
void Particles::serialize( Archive &ar, unsigned int ) {
	using boost::serialization::make_array ;

	ar & m_count ;
	ar & make_array( m_volumes.data(),      m_count ) ;
	ar & make_array( m_centers.data(),    3*m_count ) ;
	ar & make_array( m_velocities.data(), 3*m_count ) ;
	ar & make_array( m_frames.data(),     6*m_count ) ;
	ar & make_array( m_orient.data(),     6*m_count ) ;
}

}//d6

#endif

