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
	ar & make_array( m_volumes.data(),       m_count ) ;
	ar & make_array( m_centers.data(),    WD*m_count ) ;
	ar & make_array( m_velocities.data(), WD*m_count ) ;
	ar & make_array( m_frames.data(),     SD*m_count ) ;
	ar & make_array( m_orient.data(),     SD*m_count ) ;
}

template < typename Archive >
void Particles::Event::serialize( Archive &ar, unsigned int ) {
	ar & type ;
	ar & first ;
	ar & second ;
	ar & dx ;
}

}//d6

#endif

