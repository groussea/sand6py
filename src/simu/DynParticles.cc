#include "DynParticles.hh"

#include "Config.hh"
#include "Phase.hh"

#include "Scenario.hh"

#include "geo/Tensor.hh"

#include "utils/Log.hh"

#include <bogus/Core/Utils/Timer.hpp>

#include <Eigen/Eigenvalues>

#define SPLIT
#define MERGE


namespace d6 {

DynParticles::DynParticles()
{
	resize( Particles::s_MAX ) ;
}

void DynParticles::generate(const Config &c, const MeshType &mesh)
{
	m_geo.generate( Scenario::parse( c )->generator(), c.nSamples, mesh );

	 m_affine.leftCols( count() ).setZero() ;
	m_inertia.leftCols( count() ).setZero() ;

	m_meanVolume = m_geo.volumes().segment( 0, m_geo.count() ).sum() / m_geo.count() ;
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
	splitMerge( mesh );

#pragma omp parallel for
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

			grad = Du ;
//			grad.setZero() ; //FIXME
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

struct MergeInfo {
	size_t  pid ;
	Scalar  len ;
	Vec		dir ;
};

void DynParticles::splitMerge( const MeshType & mesh )
{

#ifdef SPLIT
	const std::size_t n = count() ;
	const Scalar defLength = std::pow( m_meanVolume, 1./3 ) ;

#ifdef MERGE
	std::vector< std::vector< MergeInfo > > mg_hash( mesh.nCells() ) ;
#endif

#pragma omp parallel for
	for(size_t i = 0 ; i < n ; ++i)
	{
		if( m_geo.volumes()[i] < m_meanVolume / 64 )
			continue ;
		if( m_geo.m_count + 1 == Particles::s_MAX )
			continue ;

		Mat frame ;
		tensor_view( m_geo.m_frames.col(i) ).get( frame ) ;

		Eigen::SelfAdjointEigenSolver<Mat> es(frame);
		Vec ev = es.eigenvalues().array().max(0).sqrt() ;

		typename Mat::Index kMax = 0, kMin = 0 ;
		const Scalar evMin = ev.minCoeff(&kMin) ;
		const Scalar evMax = ev.maxCoeff(&kMax) ;

		if( //ev.minCoeff() < 1.e-6
				   evMax > evMin * 4.
				&& evMax > 2 * defLength
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
//					ev[kMin] = std::max( evMin, evMax/4 ) ;

					m_geo.m_centers.col(j) = m_geo.m_centers.col(i) - ev[kMax] * es.eigenvectors().col(kMax).normalized() ;
					m_geo.m_centers.col(i) = m_geo.m_centers.col(i) + ev[kMax] * es.eigenvectors().col(kMax).normalized() ;

					m_geo.m_velocities.col(j) = m_geo.m_velocities.col(i) ;
					m_geo.m_orient.col(j) = m_geo.m_orient.col(i) ;
					m_affine.col(j) = m_affine.col(i) ;
					m_inertia(j) = m_inertia(i) ;

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
			frame = es.eigenvectors() * ev.asDiagonal() * ev.asDiagonal() * es.eigenvectors().transpose() ;
			tensor_view( m_geo.m_frames.col(i) ).set( frame ) ;

			if( evMax > 2*evMin ) {
	#ifdef MERGE
				typename MeshType::Location loc ;
				mesh.locate( m_geo.centers().col(i), loc );
				MergeInfo mgi { i, evMin, es.eigenvectors().col(kMin) }  ;
	#pragma omp critical
				mg_hash[ mesh.cellIndex( loc.cell ) ].push_back( mgi ) ;
	#endif
			}

		}

	}

	if( m_geo.m_count > Particles::s_MAX )
		m_geo.m_count = Particles::s_MAX  ;

	Log::Debug() << arg( "Split: added %1 particles, tot %2", count()-n, count() ) << std::endl ;

#ifdef MERGE

	const size_t n_before_merge = count() ;
	const size_t None = -1 ;
	std::vector< size_t > mg_indices( count(), None ) ;

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

//			m_velocities.col(i) = ( m_masses[i] * m_velocities.col(i) + m_masses[j] * m_velocities.col(j) ) / ( m_masses[i] + m_masses[j ] ) ;
//			m_affine.col(i) = ( m_masses[i] * m_affine.col(i) + m_masses[j] * m_affine.col(j) ) / ( m_masses[i] + m_masses[j ] ) ;
//			m_inertia[i] = ( m_masses[i] * m_inertia[i] + m_masses[j] * m_inertia[j] ) / ( m_masses[i] + m_masses[j ] ) ;

			m_geo.m_volumes[i] += vj ;
			m_geo.m_centers.col(i) = bary;
		}
	}


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

#endif
}

void DynParticles::resize( size_t n )
{
	m_affine.resize( 9, n );
	m_inertia.resize( 1, n );
}

} //d6
