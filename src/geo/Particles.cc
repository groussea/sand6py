#include "Particles.hh"

namespace d6 {

const size_t Particles::s_MAX = 1.e6 ;

Particles::Particles()
	: m_count(0)
{
	resize(s_MAX) ;
}

void Particles::generate(const Config &c, const MeshType &mesh)
{
	m_velocities.leftCols( count() ).setZero() ;
}

void Particles::resize(size_t n)
{
	m_volumes.resize( n );

	m_centers.resize( 3, n);
	m_velocities.resize( 3, n);

	m_frames.resize( 6, n);
	m_orient.resize( 6, n);
}

} //d6
