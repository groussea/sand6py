#ifndef D6_PARTICLES_IO_HH
#define D6_PARTICLES_IO_HH

#include "Particles.hh"

#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>

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

	ar & m_log ;
}

template < typename Archive >
void Particles::Event::serialize( Archive &ar, unsigned int ) {
	ar & type ;
	ar & first ;
	ar & second ;
}

}//d6

#endif

