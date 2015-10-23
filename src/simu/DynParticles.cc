#include "DynParticles.hh"

#include "utils/Config.hh"

namespace d6 {

DynParticles::DynParticles()
{
	resize( Particles::s_MAX ) ;
}

void DynParticles::generate(const Config &c, const MeshType &mesh)
{
	m_geo.generate( c, mesh );

	m_geo.m_masses *= c.volMass ;

	m_affine.leftCols( count() ).setZero() ;
}

void DynParticles::resize( size_t n )
{
	m_affine.resize( 9, n );
}

} //d6
