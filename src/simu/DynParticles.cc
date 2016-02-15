#include "DynParticles.hh"

#include "Phase.hh"

#include "Scenario.hh"

#include "geo/Tensor.hh"

#include "utils/Log.hh"
#include "utils/Config.hh"

#include <bogus/Core/Utils/Timer.hpp>

#include <Eigen/Eigenvalues>

#include <random>

//TODO allow resampling particles w/ unstructured dual
#ifndef D6_UNSTRUCTURED_DUAL
#define SPLIT
#define MERGE
#endif

#define GRAD_FROM_VEL

namespace d6 {

DynParticles::DynParticles()
{
	resize( Particles::s_MAX ) ;
}

void DynParticles::generate(const Config &c, const MeshType &mesh, const Scenario &scenario)
{
	//Generate particles according to grid

	m_geo.generate( scenario.generator(), c.nSamples, mesh, c.cohesion > 0, c.initialOri );

	  m_affine.leftCols( count() ).setZero() ;
	 m_inertia.leftCols( count() ).setZero() ;
	m_cohesion.leftCols( count() ).setConstant( 1. ) ;

	m_meanVolume = m_geo.volumes().segment( 0, m_geo.count() ).sum() / m_geo.count() ;

	// Randomize particle positions
	if( c.randomize > 0 ) {

		std::uniform_real_distribution<Scalar> dist( -c.randomize, c.randomize );

#pragma omp parallel
		{
			std::random_device rd;
			std::default_random_engine e( rd() ) ;

#pragma omp for
			for( size_t i = 0 ; i < count() ; ++i ) {
				const Scalar size = std::sqrt( m_geo.frames().col(i)[0] ) ;
				Vec3 rd_depl( dist(e), dist(e), dist(e) ) ;
				m_geo.m_centers.col(i) += size * rd_depl.head<WD>() ;
				clamp_particle( i, mesh );
			}

		}
	}

}

void DynParticles::clamp_particle(size_t i, const MeshType &mesh)
{
	const Vec p = mesh.clamp_point( m_geo.centers().col(i) ) ;
	const Vec dp = m_geo.centers().col(i) - p ;
	const Scalar dpn2 = p.squaredNorm() ;

	if( dpn2 > 1.e-12 ) {
		// Particle was outside domain, project-out normal velocity
		m_geo.m_velocities.col(i) -= ( dp.dot(m_geo.m_velocities.col(i)) ) * dp / dpn2 ;
	}

	m_geo.m_centers.col(i) = p ;
}

void DynParticles::update(const Config &config, const Scalar dt, const Phase &phase )
{
	bogus::Timer timer ;

	const std::size_t n = count() ;

	const PrimalShape& shape = phase.velocity.shape() ;
	const MeshType& mesh = shape.derived().mesh() ;

	// Split and merge particles before setting their velocities / moving them
	splitMerge( mesh );

#pragma omp parallel for
	for( size_t i = 0 ; i < n ; ++i ) {

		const Vec &p0 =  m_geo.m_centers.col(i) ;

		typename PrimalShape::Location p0loc ;
		mesh.locate( p0, p0loc );

		// Grid velocity at current position
		const Vec  vi ( phase.velocity(p0loc) ) ;

		m_geo.m_velocities.col(i) = vi ; //PIC
		m_geo.m_centers.col(i) += phase.geo_proj(p0loc) + dt * ( m_geo.m_velocities.col(i) ) ;

		clamp_particle( i, mesh );

		//APIC
		Mat grad ;
		{
			// Recompute gradient from velocities to avoid smoothing

			typename PrimalShape::Interpolation itp ;
			typename PrimalShape::Derivatives derivatives ;

			shape.interpolate( p0loc, itp );
			shape.get_derivatives( p0loc, derivatives ) ;

			Mat Bp = Mat::Zero() ;

			for( Index k = 0 ; k < itp.nodes.rows() ; ++k ) {
				Bp += phase.velocity[ itp.nodes[k] ] * derivatives.row(k) ;
			}

			tensor_view( m_affine.col(i) ).set( Bp );

			// Debonding function
			m_cohesion(i) /= ( 1. + dt * config.cohesion_decay * (Bp+Bp.transpose()).norm() ) ;

#ifdef GRAD_FROM_VEL
			grad = Bp ;
#endif
		}


		// Dual Quantities
		typename DualShape::Location d0loc ;
		phase.stresses.shape().locate_by_pos_or_id( p0, i, d0loc ) ;

		// Frames and orientation
		Mat Du, Wu ;
		phase.sym_grad.get_sym_tensor( d0loc, Du );
		phase.spi_grad.get_spi_tensor( d0loc, Wu );

		// Inertia
		{
			const Scalar DuT = ( Du - 1./WD * Du.trace() * Mat::Identity() ).norm()  ;
			m_inertia(i) = DuT / std::sqrt( std::max( 1.e-16, phase.stresses(d0loc)[0] ) ) ;
		}

		// Frame
		{
#ifndef GRAD_FROM_VEL
			const Mat grad = Wu + Du ;
#endif

			auto  frame_view( tensor_view( m_geo.m_frames.col(i) ) ) ;
			Mat frame ;
			frame_view.get( frame ) ;
			frame += dt * ( grad * frame + frame * grad.transpose() ) ;
			frame_view.set( frame ) ;
		}

		//Orientation
		{
			const Scalar lambda = config.elongation ;

			auto orient_view( tensor_view( m_geo.m_orient.col(i) ) ) ;
			Mat orient ;
			orient_view.get( orient ) ;

			Eigen::SelfAdjointEigenSolver<Mat> es(orient);
			Vec ev = es.eigenvalues().array().max(0).min(1) ;

			// Translate to tangent space
			Mat a2 = .5 * es.eigenvectors() * (1. - ev.array()).matrix().asDiagonal() * es.eigenvectors().transpose() ;

			// Quadratic closure
			Mat a4grad = (a2.cwiseProduct( Du ).sum() * a2 ) ;

			a2 += dt * ( Wu * a2  - a2 * Wu + lambda * ( Du * a2  + a2 * Du ) ) ;
			a2 -= 2 * lambda * dt  * a4grad ;

			//Back to normal
			es.compute(a2);
			ev = es.eigenvalues().array().max(0).min(1) ;

			orient = es.eigenvectors() * (1 - (WD-1)*ev.array()).matrix().asDiagonal() * es.eigenvectors().transpose() ;

			// Folgar-Tucker brownian motion
			orient += dt * config.brownian * Du.norm() * ( Mat::Identity() - WD * orient ) ;

			// Normalize
			orient /= orient.lpNorm<1>() ;

			orient_view.set( orient ) ;
		}

	}

	Log::Debug() << "Particles advection time: " << timer.elapsed() << std::endl ;

}

void DynParticles::integratePrimal( std::vector< bool > &activeCells,
			   PrimalScalarField &phi,    PrimalVectorField &phiVel) const
{
	typedef typename PrimalScalarField::ShapeFuncImpl Shape ;

	const Shape& shape = phi.shape().derived() ;
	activeCells.assign( shape.mesh().nCells(), false );

	phi.set_zero();
	phiVel.set_zero();

	for ( size_t i = 0 ; i < count() ; ++i ) {

		const Scalar m = m_geo.volumes()[i] ;
		const Vec p0 = m_geo.centers().col(i) ;

		typename Shape::Location loc ;
		shape.locate( p0, loc );

		typename Shape::Interpolation itp ;
		shape.interpolate_tpz( loc, itp );

		activeCells[ shape.mesh().cellIndex( loc.cell ) ] = true ;

		// Accumulate particle data at grid nodes
		phi.add_at( itp, m );
		phiVel .add_at( itp, m * m_geo.velocities().col(i) );


		// APIC
		{
			Mat affine ;
			tensor_view( m_affine.col(i) ).get( affine ) ;

			typename MeshType::CellGeo cellGeo ;
			shape.mesh().get_geo( loc.cell, cellGeo );

			for( Index k = 0 ; k < MeshType::NV ; ++k ) {
				phiVel[ itp.nodes[k] ] += m * itp.coeffs[k] * ( affine * ( cellGeo.vertex( k ) - p0 ) ) ;
			}
		}
	}
}

void DynParticles::integrateDual(
								  DualScalarField &phi,      DualScalarField &phiInertia,
								  DualTensorField &phiOrient,DualScalarField &phiCohesion
								  ) const
{
	typedef DualShape Shape ;
	const Shape& shape = phi.shape().derived() ;

	phi.set_zero();
	phiInertia .set_zero();
	phiCohesion.set_zero();
	phiOrient  .set_zero();

	for ( size_t i = 0 ; i < count() ; ++i ) {

		const Scalar m = m_geo.volumes()[i] ;
		const Vec p0 = m_geo.centers().col(i) ;

		typename Shape::Location loc ;
		shape.locate_by_pos_or_id( p0, i, loc );

		typename Shape::Interpolation itp ;
		shape.interpolate( loc, itp );

		phi.add_at( itp, m );
		phiInertia .add_at( itp, m * m_inertia[i] );
		phiCohesion.add_at( itp, m * m_cohesion[i] );
		phiOrient  .add_at( itp, m * m_geo.orient().col(i) );
	}

}

struct MergeInfo {
	size_t  pid ;
	Scalar  len ;
	Vec		dir ;
};

void DynParticles::splitMerge( const MeshType & mesh )
{

#ifdef SPLIT
	const std::size_t n = count() ;
	const Scalar defLength = std::pow( m_meanVolume, 1./WD ) ;

#ifdef MERGE
	std::vector< std::vector< MergeInfo > > mg_hash( mesh.nCells() ) ;
#endif

#pragma omp parallel for
	for(size_t i = 0 ; i < n ; ++i)
	{

		Mat frame ;
		tensor_view( m_geo.m_frames.col(i) ).get( frame ) ;

		Eigen::SelfAdjointEigenSolver<Mat> es(frame);
		Vec ev = es.eigenvalues().array().max(0).sqrt() ;

		typename Mat::Index kMax = 0, kMin = 0 ;
		const Scalar evMin = ev.minCoeff(&kMin) ;
		const Scalar evMax = ev.maxCoeff(&kMax) ;

		if(  	   evMax > evMin * 4.     // Eigenvalues ratio
				&& evMax > defLength      // Avoid splitting too small particles
				&& m_geo.volumes()[i] > m_meanVolume / 64 // Avoid splitting too ligth particles
				&& m_geo.m_count + 1 != Particles::s_MAX // Avoid running out of memory
				   )
		{
			// Split
			size_t j = 0 ;
#pragma omp atomic capture
			j = m_geo.m_count++ ;

			if( j < Particles::s_MAX ) {

					m_geo.m_volumes[i] *= .5 ;
					m_geo.m_volumes[j] = m_geo.m_volumes[i] ;

					ev[kMax] *= .5 ;
					ev[kMin] = std::max( defLength / 8, ev[kMin] ) ;

					m_geo.m_centers.col(j) = m_geo.m_centers.col(i) - ev[kMax] * es.eigenvectors().col(kMax).normalized() ;
					m_geo.m_centers.col(i) = m_geo.m_centers.col(i) + ev[kMax] * es.eigenvectors().col(kMax).normalized() ;

					m_geo.m_velocities.col(j) = m_geo.m_velocities.col(i) ;
					m_geo.m_orient.col(j) = m_geo.m_orient.col(i) ;
					m_affine.col(j) = m_affine.col(i) ;
					m_inertia(j) = m_inertia(i) ;

					m_cohesion(j) = m_cohesion(i) ;

					clamp_particle( i, mesh );
					clamp_particle( j, mesh );

					frame = es.eigenvectors() * ev.asDiagonal() * ev.asDiagonal() * es.eigenvectors().transpose() ;

					tensor_view( m_geo.m_frames.col(i) ).set( frame ) ;
					tensor_view( m_geo.m_frames.col(j) ).set( frame ) ;

					const Vec dx = .5 * (m_geo.m_centers.col(j) - m_geo.m_centers.col(i)) ;
					m_events.log( Particles::Event::split( i, j, dx ) );

			}
		} else {

			//Repair flat frames
			ev[kMin] = std::max( defLength / 8, ev[kMin] ) ;
			ev[kMax] = std::min( defLength * 8, ev[kMax] ) ;
			frame = es.eigenvectors() * ev.asDiagonal() * ev.asDiagonal() * es.eigenvectors().transpose() ;
			tensor_view( m_geo.m_frames.col(i) ).set( frame ) ;

			// Add particle to merge candidates
	#ifdef MERGE
			if( evMax > 2*evMin ) {
				typename MeshType::Location loc ;
				mesh.locate( m_geo.centers().col(i), loc );
				MergeInfo mgi { i, evMin, es.eigenvectors().col(kMin) }  ;
	#pragma omp critical
				mg_hash[ mesh.cellIndex( loc.cell ) ].push_back( mgi ) ;
			}
	#endif

		}

	}

	if( m_geo.m_count > Particles::s_MAX )
		m_geo.m_count = Particles::s_MAX  ;

	Log::Debug() << arg( "Split: added %1 particles, tot %2", count()-n, count() ) << std::endl ;

#ifdef MERGE

	const size_t n_before_merge = count() ;
	const size_t None = -1 ;
	std::vector< size_t > mg_indices( count(), None ) ;

	// Look for overlapping particles in each mesh cell
#pragma omp parallel for
	for( size_t cidx = 0 ; cidx < mg_hash.size() ; ++ cidx ) {
		const std::vector< MergeInfo >& list = mg_hash[cidx] ;
		const unsigned m = list.size() ;
		for( unsigned k = 0 ; k < m ; ++k ) {
			if( mg_indices[list[k].pid] != None ) continue ;
			for( unsigned l = 0 ; l < k ; ++ l ) {
				if( mg_indices[list[l].pid] != None ) continue ;

				const Vec pk = m_geo.centers().col( list[k].pid ) ;
				const Vec pl = m_geo.centers().col( list[l].pid ) ;

				const Scalar depl = std::fabs( list[k].dir.dot( list[l].dir ) ) * ( list[k].len + list[l].len ) ;

				if( (pk - pl).squaredNorm() < depl*depl ) {
					mg_indices[ list[k].pid ] = list[l].pid ;
					mg_indices[ list[l].pid ] = list[k].pid ;

					break ;
				}
			}
		}
	}

	// Compute result of merges
#pragma omp parallel for
	for( size_t i = 0 ; i < mg_indices.size() ; ++ i ) {
		if ( mg_indices[i] != None && mg_indices[i]>i ) {
			const size_t j = mg_indices[i] ;

			const Vec pi = m_geo.centers().col( i ) ;
			const Vec pj = m_geo.centers().col( j ) ;
			const Scalar vi = m_geo.volumes()[i] ;
			const Scalar vj = m_geo.volumes()[j] ;
			const Vec bary = ( vi * pi + vj * pj ) / ( vi + vj ) ;

			Mat fi, fj ;
			tensor_view( m_geo.m_frames.col(i) ).get( fi ) ;
			tensor_view( m_geo.m_frames.col(j) ).get( fj ) ;

			Mat frame = ( ( fi + (pi - bary)*(pi - bary).transpose() ) * vi +
						  ( fj + (pj - bary)*(pj - bary).transpose() ) * vj )
					/ ( vi + vj ) ;
			tensor_view( m_geo.m_frames.col(i) ).set( frame ) ;

			m_geo.m_orient.col(i) = ( vi * m_geo.orient().col(i) + vj * m_geo.orient().col(j) ) / ( vi + vj ) ;

			m_cohesion[i] = ( vi * m_cohesion[i] + vj * m_cohesion[j] ) / ( vi + vj ) ;

			m_geo.m_volumes[i] += vj ;
			m_geo.m_centers.col(i) = bary;
		}
	}

	// Reallocate particles (ensure contiguous indices)
	for( size_t j = 0 ; j < n_before_merge ; ++ j ) {
		if ( mg_indices[j] != None && mg_indices[j] < j ) {
			//i is empty

			const Vec dx = m_geo.m_centers.col( mg_indices[j] ) - m_geo.m_centers.col( j ) ;
			m_events.log( Particles::Event::merge( mg_indices[j], j, dx ) );

			if( j >= m_geo.m_count ) continue ;

			size_t reloc_src = -1 ;

			// Find last non-emptym pos
			do {
				reloc_src = --m_geo.m_count ;
			} while( reloc_src > j && mg_indices[reloc_src] != None
					 && mg_indices[reloc_src] < reloc_src  ) ;

			m_geo.m_volumes[j] = m_geo.m_volumes[reloc_src] ;
			m_geo.m_centers.col(j) = m_geo.m_centers.col(reloc_src) ;
			m_geo.m_orient.col(j) = m_geo.m_orient.col(reloc_src) ;
			m_geo.m_frames.col(j) = m_geo.m_frames.col(reloc_src) ;

		}

	}


	Log::Debug() << arg( "Merge: removed %1 particles, tot %2", (n_before_merge-count()), count() ) << std::endl ;

#endif

#else
	(void) mesh ;
#endif
}

void DynParticles::remove(size_t j)
{
	size_t src = 0 ;
#pragma omp critical
{
	src = --m_geo.m_count ;
	m_events.log( Particles::Event::remove(j, src) );

	if( j < src ) {
		m_geo.m_volumes[j] = m_geo.m_volumes[src] ;
		m_geo.m_centers.col(j) = m_geo.m_centers.col(src) ;
		m_geo.m_velocities.col(j) = m_geo.m_velocities.col(src) ;
		m_geo.m_orient.col(j) = m_geo.m_orient.col(src) ;
		m_geo.m_frames.col(j) = m_geo.m_frames.col(src) ;

		m_inertia.col(j) = m_inertia.col(src) ;
		m_cohesion.col(j) = m_cohesion.col(src) ;
		m_affine.col(j) = m_affine.col(src) ;
	}
}

}

void DynParticles::resize( size_t n )
{
	m_affine.resize( WD*WD, n );
	m_inertia.resize( 1, n );
	m_cohesion.resize( 1, n );
}

} //d6
