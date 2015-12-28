#include "Sampler.hh"

#include "Offline.hh"

#include "geo/Tensor.hh"
#include "simu/Phase.hh"

#include <Eigen/Eigenvalues>

#include <random>
#include <iostream>

#define RAD_FAC 3

namespace d6 {

struct BallSampler {

	BallSampler()
	: e(rd()), dist(-1,1)
	{}

	void operator() (Vec &v) {
		do {
			v << dist(e), dist(e), dist(e) ;
		} while (v.squaredNorm() > 1. ) ;
	}

	void onSphere( Vec& v ) {
		do {
			(*this)(v) ;
		} while( v.isZero() ) ;

		v /= v.norm() ;
	}

	std::random_device rd;
	std::default_random_engine e ;
	std::uniform_real_distribution<Scalar> dist;
} ;

void Sampler::move( const Scalar dt )
{
	// offset += DU*offset * st

	const Phase& grains = m_offline.grains() ;

	const Index n = count() ;

#pragma omp parallel for
	for(Index i = 0 ; i < n ; ++i) {

		const Vec p0 = m_positions.col(i).cast< Scalar >()  ;
		const Vec p1 = p0 + dt * grains.velocity( p0 ) ;
		m_positions.col(i) = p1.cast< float >() ;
	}

	// Predict particle positions
	const Index np = m_particlesCount ;
#pragma omp parallel for
	for(Index i = 0 ; i < np ; ++i) {
		m_predPos.col(i) += dt*grains.velocity(	m_predPos.col(i) )   ;
	}

}

void Sampler::compute_absolute( )
{
	// offset += DU*offset * st


	const Scalar noise = 1. ;

	const Particles& particles = m_offline.particles() ;
	const Phase& grains = m_offline.grains() ;

	assert( particles.count() == m_particlesCount ) ;

	const Index n = count() ;

#pragma omp parallel for
	for( Index i = 0 ; i < n ; ++i) {
		const unsigned pid = m_particleIds[i] ;

		const Vec p0 = particles.centers().col( pid ) ;
		const Vec p0_pred = m_predPos.col( pid ) ;

		m_offsets.col( i ) = m_positions.col(i).cast<Scalar>() - p0_pred ;

		Mat frame ;
		tensor_view( particles.frames().col(pid) ).get( frame ) ;
		frame *= RAD_FAC ;

		const Scalar vn = std::max( 1., Scalar( m_offsets.col(i).transpose() *  frame.inverse() * m_offsets.col(i) ) ) ;
		m_offsets.col( i ) /= std::sqrt(vn)  ;


		const Vec pos = (p0 + m_offsets.col(i) ) ;
		m_positions.col(i) = pos.cast< float >() ;

		Eigen::Vector3f grad_phi = grains.grad_phi( pos ).cast<float>() ;
		grad_phi += 1.e-2 * Eigen::Vector3f( 0, 0, -1 ) ;
		const float gn = grad_phi.norm() ;
		if( gn > 1.e-6 )
			grad_phi /= gn ;
		else
			grad_phi = Eigen::Vector3f( 0, 0, -1 ) ;


		Eigen::Quaternionf rot( Eigen::AngleAxisf( noise, m_normalNoise.col(i) ) ) ;
		m_normals.col(i) =  - ( rot * grad_phi ).normalized() ;

		if( m_visibility(i) >= 0 ) {
			m_visibility(i) = std::max( 0., std::min( 1., 1.5 - 1.5*grains.fraction(pos) ) ) ;
		}
	}

}

void Sampler::reassign( )
{
	typedef Particles::Event Event ;

	const unsigned NoEvent = -1 ;
	const unsigned Destroy = -2 ;

	std::vector< unsigned > eventIds ;

	const Particles& particles = m_offline.particles() ;

	// Reset particle positions and reserve storage space
	unsigned nEvents = 0 ;
	for( const std::vector< Event > &events :  m_offline.events().log() ) {
		nEvents += events.size() ;
	}
	m_predPos.resize(3, m_particlesCount + nEvents );
	m_predPos.leftCols( m_particlesCount ) = particles.centers().leftCols( m_particlesCount ) ;

	const Index n = count() ;
	unsigned nParticles = m_particlesCount ;

	for(const std::vector< Event > &events :  m_offline.events().log()) {

		unsigned splits = 0 ;
		unsigned merges = 0 ;
		unsigned removes = 0 ;

		// Assign event-ids to particles
		eventIds.assign( nParticles + events.size(), NoEvent );
		for( unsigned eId = 0 ; eId < events.size() ; ++eId ) {
			const Event& e = events[eId] ;

			if( e.type == Event::Split ) {
				++splits ;

				eventIds[ e.first  ] = eId ;

				m_predPos.col( e.second ) = m_predPos.col( e.first ) + e.dx ;
				m_predPos.col( e.first  ) = m_predPos.col( e.first ) - e.dx ;
			} else if( e.type == Event::Remove ) {
				++removes ;
			} else if( e.type == Event::Merge ) {
				++merges;

				eventIds[ e.first  ] = eId ;
				eventIds[ e.second ] = eId ;

				m_predPos.col( e.first  ) = m_predPos.col( e.second ) + e.dx ;
			}
		}

		//Split -- merge
#pragma omp parallel for
		for( Index i = 0 ; i < n ; ++i ) {

			const unsigned pid = m_particleIds[i] ;
			if( eventIds[ pid ] == NoEvent ) continue ;

			const Event& e = events[eventIds[pid]] ;

			if( e.type == Event::Split ) {

				const Vec dx = events[eventIds[pid]].dx ;

				//dx = new_center_2 - old_center = old_center - new_center_1
				//offset = p - old_center
				//offset_2 = p - new_center_2 = p - (dx + old_center) = p - dx
				//offset_1 = p - new_center_1 = p - (old_center - dx) = p + dx

				const Vec &off1 = m_offsets.col(i) + dx ;
				const Vec &off2 = m_offsets.col(i) - dx ;

				if( off1.squaredNorm() < off2.squaredNorm() ) {
					m_offsets.col(i) = off1 ;
				} else {
					m_offsets.col(i) = off2 ;
					m_particleIds[i] = e.second ;
				}
			}
		}

#pragma omp parallel for
		for( Index i = 0 ; i < n ; ++i ) {

			const unsigned pid = m_particleIds[i] ;
			if( eventIds[ pid ] == NoEvent ) continue ;

			const Event& e = events[eventIds[pid]] ;

			if( e.type == Event::Merge ) {

				const Vec dx = events[eventIds[pid]].dx ;
				if( pid == e.first ) {
					m_offsets.col(i) -= dx ;
				} else {
					m_offsets.col(i) += dx ;
					m_particleIds[i] = e.first ;
				}
			}
		}

		nParticles += splits ;

		// Reloc
		std::vector< unsigned > relocs ;
		relocs.assign( nParticles, NoEvent );
		size_t reloc_src = nParticles ;
		for( unsigned eId = 0 ; eId < events.size() ; ++eId ) {
			const Event& e = events[eId] ;

			if( e.type == Event::Merge && e.second < nParticles ) {
				// Find last non-emptym pos
				do {
					--reloc_src ;
				} while( reloc_src > e.second && eventIds[reloc_src] != NoEvent
						 && events[ eventIds[reloc_src] ].type == Event::Merge
						 && events[ eventIds[reloc_src] ].second == reloc_src
						 ) ;

				relocs[ reloc_src ] = e.second ;
				m_predPos.col( e.second ) = m_predPos.col( reloc_src ) ;
			}
		}

#pragma omp parallel for
		for( Index i = 0 ; i < n ; ++i ) {

			const unsigned pid = m_particleIds[i] ;
			if( relocs[ pid ] == NoEvent ) continue ;

			m_particleIds[ i ] = relocs[pid] ;
		}

		nParticles -= merges ;

		// Remove -- Impossible without destroyign samples
		relocs.assign( nParticles, NoEvent );
		for( unsigned eId = 0 ; eId < events.size() ; ++eId ) {
			const Event& e = events[eId] ;

			if( e.type == Event::Remove ) {
				relocs[ e.second ] = e.first ;
				relocs[ e.first ] = Destroy ;
				m_predPos.col( e.first ) = m_predPos.col( e.second ) ;
			}
		}
#pragma omp parallel for
		for( Index i = 0 ; i < n ; ++i ) {

			const unsigned pid = m_particleIds[i] ;
			if( relocs[ pid ] == NoEvent ) continue ;

			if( relocs[ pid ] == Destroy )  {
				m_particleIds[ i ] = 0 ;
				m_visibility(i) = -1 ;
				continue ;
			}

			m_particleIds[ i ] = relocs[pid] ;
		}

		nParticles -= removes ;

	}

	m_particlesCount = nParticles ;

}

void Sampler::sampleParticles( unsigned nSamples )
{
	const Particles& particles = m_offline.particles() ;
	const Index n = particles.count() ;

	m_particleIds.resize( n * nSamples ) ;
	m_offsets.resize( 3, n * nSamples ) ;
	m_normalNoise.resize( 3, n * nSamples ) ;

	m_positions.resize( 3, n * nSamples ) ;
	m_normals.resize( 3, n * nSamples ) ;

	m_visibility.resize( n * nSamples ) ;
	m_visibility.setZero() ;

#pragma omp parallel
	{
		BallSampler sampler ;

		#pragma omp for
		for(Index i = 0 ; i < n ; ++i) {

			Vec v ;
			Mat frame ;
			tensor_view( particles.frames().col(i) ).get( frame ) ;
			frame *= RAD_FAC ;

			for( unsigned j = 0 ; j < nSamples ; ++j )
			{
				const unsigned idx = i*nSamples +  j ;
				m_particleIds[idx] = i ;

				sampler(v) ;

				const Scalar vn = std::max( 1., Scalar( v.transpose() * frame.inverse() * v ) ) ;
				m_offsets.col( idx ) = v/std::sqrt(vn)  ;
				m_positions.col( idx ) = ( particles.centers().col(i) + m_offsets.col( idx )).cast<float>() ;

				sampler( v ) ;

				m_normalNoise.col( idx ) = v.cast< float >() ;

			}

		}
	}

	m_particlesCount = n ;
	m_predPos = particles.centers().leftCols( m_particlesCount ) ;

	compute_absolute() ;
}


} // d6
