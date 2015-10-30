#include "DynParticles.hh"

#include "Config.hh"
#include "Phase.hh"

#include "geo/Tensor.hh"

#include "utils/Log.hh"
#include <bogus/Core/Utils/Timer.hpp>

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

void DynParticles::clamp_particle(size_t i, const MeshType &mesh)
{
	const Vec p = mesh.clamp_point( m_geo.centers().col(i) ) ;
	const Vec dp = m_geo.centers().col(i) - p ;
	const Scalar dpn2 = p.squaredNorm() ;

	if( dpn2 > 1.e-12 ) {
		// Non-zero projection, project-out normal velocity
		m_geo.m_velocities.col(i) -= ( dp.dot(m_geo.m_velocities.col(i)) ) * dp / dpn2 ;
	}

	m_geo.m_centers.col(i) = p ;
}

void DynParticles::update(const Config &config, const Phase &phase )
{
	bogus::Timer timer ;

	const Scalar dt = config.dt() ;
	const std::size_t n = count() ;

	const MeshType& mesh = phase.velocity.mesh() ;

	for( size_t i = 0 ; i < n ; ++i ) {

		const Vec &p0 =  m_geo.m_centers.col(i) ;

		const Vec  vi ( phase.velocity(p0) ) ;

		m_geo.m_velocities.col(i) = vi ; //PIC
		m_geo.m_centers.col(i) += phase.geo_proj(p0) + dt * ( m_geo.m_velocities.col(i) ) ;

		clamp_particle( i, mesh );

		//APIC
		{
			// Recompute gradient from velocities to avoid smoothing

			typename MeshType::Location loc ;
			typename MeshType::Interpolation itp ;
			typename MeshType::Derivatives derivatives ;

			mesh.locate( p0, loc );
			mesh.interpolate( loc, itp );
			mesh.get_derivatives( loc, derivatives ) ;

			Mat Bp = Mat::Zero() ;

			for( Index k = 0 ; k < MeshType::NV ; ++k ) {
				Bp += phase.velocity[ itp.nodes[k] ] * derivatives.row(k) ;
			}

			tensor_view( m_affine.col(i) ).set( Bp );
		}


		// Frames and orientation
		Mat grad ;

		// Inertia
		{
			Mat Du ;
			phase.sym_grad.get_sym_tensor( p0, Du );

			//FIXME
			//grad = Du ;
			grad.setZero() ;
			phase.spi_grad.add_spi_tensor( p0, grad );

			const Scalar DuT = ( Du - 1./3. * Du.trace() * Mat::Identity() ).norm()  ;
			m_inertia(i) = DuT / std::sqrt( std::max( 1.e-16, phase.stresses(p0)[0] ) ) ;

		}

		// Frame
		{
			auto  frame_view( tensor_view( m_geo.m_frames.col(i) ) ) ;
			Mat frame ;
			frame_view.get( frame ) ;
			frame += dt * ( grad * frame + frame * grad.transpose() ) ;
			frame_view.set( frame ) ;
		}

		//Orientation
		{
			auto orient_view( tensor_view( m_geo.m_orient.col(i) ) ) ;
			Mat a2 ;
			orient_view.get( a2 ) ;

			// Quadratic closure
			Mat a4grad = (a2.cwiseProduct( grad ).sum() * a2 ) ;

			a2    += dt * ( grad * a2    +    a2 * grad.transpose() ) ;

			a2 -= 2 * dt  * a4grad ;
			a2 = a2.normalized() ;

			orient_view.set( a2 ) ;
		}

	}

	Log::Debug() << "Particles advection time: " << timer.elapsed() << std::endl ;

}

void DynParticles::read(std::vector<bool> &activeCells,
						ScalarField &phi, VectorField &phiVel,
						ScalarField &phiInertia, TensorField &phiOrient) const
{
	const MeshType& mesh = phi.mesh() ;
	activeCells.assign( mesh.nCells(), false );

	phi.set_zero();
	phiVel.set_zero();
	phiInertia.set_zero();
	phiOrient.set_zero();

	for ( size_t i = 0 ; i < count() ; ++i ) {

		const Scalar m = m_geo.volumes()[i] ;
		const Vec p0 = m_geo.centers().col(i) ;

		typename MeshType::Location loc ;
		typename MeshType::Interpolation itp ;
		mesh.locate( p0, loc );
		mesh.interpolate( loc, itp );

		activeCells[ mesh.cellIndex( loc.cell ) ] = true ;

			   phi.add_at( itp, m );
			phiVel.add_at( itp, m * m_geo.velocities().col(i) );
		phiInertia.add_at( itp, m * m_inertia[i] );
		phiOrient .add_at( itp, m * m_geo.orient().col(i) );


		// APIC
		{
			Mat affine ;
			tensor_view( m_affine.col(i) ).get( affine ) ;

			typename MeshType::CellGeo cellGeo ;
			mesh.get_geo( loc.cell, cellGeo );

			for( Index k = 0 ; k < MeshType::NV ; ++k ) {
				phiVel[ itp.nodes[k] ] += m * itp.coeffs[k] * ( affine * ( cellGeo.vertex( k ) - p0 ) ) ;
			}
		}
	}
}

void DynParticles::resize( size_t n )
{
	m_affine.resize( 9, n );
	m_inertia.resize( 1, n );
}

} //d6
