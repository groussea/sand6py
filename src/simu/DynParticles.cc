#include "DynParticles.hh"

#include "simu/Config.hh"

namespace d6 {

DynParticles::DynParticles()
{
	resize( Particles::s_MAX ) ;
}

void DynParticles::generate(const Config &c, const MeshType &mesh)
{
	auto phi = [&]( const Vec& x){ return x[2] > .5*c.box[2] ? 1. : 0. ; } ;

	m_geo.generate( make_expr( phi ), c.nSamples, mesh );

	m_affine.leftCols( count() ).setZero() ;
}

void DynParticles::resize( size_t n )
{
	m_affine.resize( 9, n );
}

} //d6
