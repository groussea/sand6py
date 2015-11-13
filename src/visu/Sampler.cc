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

	const Particles& particles = m_offline.particles() ;
	const Phase& grains = m_offline.grains() ;

	const Index n = count() ;

#pragma omp parallel for
	for(Index i = 0 ; i < n ; ++i) {
		const unsigned pid = m_particleIds[i] ;

//		const Vec p0 = particles.centers().col( pid ) ;
//		const Vec p0 = m_positions.col(i).cast< Scalar >() - m_offsets.col(i) ;

//		Mat grad ;
//		tensor_view( grains.sym_grad( p0 ) ).get( grad ) ;
//		Mat spin ;
//		tensor_view( grains.spi_grad( p0 ) ).get( spin ) ;

//		grad += spin ;
//		m_offsets.col(i) += dt * grad * m_offsets.col(i) ;

//		m_normalNoise.col(i) = ( m_normalNoise.col(i) + dt * spin * m_normalNoise.col(i) ).normalized() ;

		const Vec p0 = m_positions.col(i).cast< Scalar >()  ;
		const Vec p1 = p0 + dt * grains.velocity( p0 ) ;
		m_positions.col(i) = p1.cast< float >() ;
	}

}

void Sampler::compute_absolute( )
{
	// offset += DU*offset * st

	const Scalar noise = 1. ;

	const Particles& particles = m_offline.particles() ;
	const Phase& grains = m_offline.grains() ;

	const Index n = count() ;

#pragma omp parallel for
	for( Index i = 0 ; i < n ; ++i) {
		const unsigned pid = m_particleIds[i] ;

		const Vec p0 = particles.centers().col( pid ) ;
		m_offsets.col( i ) = m_positions.col(i).cast<Scalar>() - p0 ;

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

		m_visibility(i) = std::max( 0., std::min( 1., 1.5 - 1.5*grains.fraction(pos) ) ) ;
	}

}

void Sampler::reassign( )
{
	typedef Particles::Event Event ;

	const unsigned NoEvent = -1 ;

	std::vector< unsigned > eventIds ;

	const Particles& particles = m_offline.particles() ;

	unsigned nParticles = m_particlesCount ;

	const Index n = count() ;

	for( unsigned s = 0 ; s < m_offline.log().log().size() ; ++s ) {
		const std::vector< Event > &events = m_offline.log().log()[s] ;

		unsigned splits = 0 ;
		unsigned merges = 0 ;

		// Assign event-ids to particles
		eventIds.assign( n, NoEvent );
		for( unsigned eId = 0 ; eId < events.size() ; ++eId ) {
			const Event& e = events[eId] ;
			if( e.type == Event::Split ) {
				eventIds[ e.first  ] = eId ;
				++splits ;
			}
		}

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

		nParticles += splits - merges ;

	}

	m_particlesCount = nParticles ;

	std::cout << "Npartic;es == " << m_particlesCount << std::endl ;
	std::cout << "Nsampeses == " << count() << std::endl ;

	// partcles Ids <-> SM
	// offsets <->

	// TOD0
	// Require 3 pos
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
	compute_absolute() ;
}


} // d6
