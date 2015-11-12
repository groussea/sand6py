#include "Sampler.hh"

#include "Offline.hh"

#include "geo/Tensor.hh"
#include "simu/Phase.hh"

#include <Eigen/Eigenvalues>
#include <random>

namespace d6 {

struct BallSampler {

	BallSampler()
	: e(rd()), dist(0,1)
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

void Sampler::updateOffsets( const Scalar dt )
{
	// offset += DU*offset * st

	const Particles& particles = m_offline.particles() ;
	const Phase& grains = m_offline.grains() ;

	const Index n = count() ;
	
#pragma omp for
	for(Index i = 0 ; i < n ; ++i) {
		const unsigned pid = m_particleIds[i] ;
		const Vec p0 = particles.centers().col( pid ) ;

		Mat grad ; 
		tensor_view( grains.sym_grad( p0 ) ).get( grad ) ;
		Mat spin ; 
		tensor_view( grains.spi_grad( p0 ) ).get( spin ) ;

		grad += spin ;
		m_offsets.col(i) += dt * grad * m_offsets.col(i) ;
		
		m_normalNoise.col(i) = ( m_normalNoise.col(i) + dt * spin * m_normalNoise.col(i) ).normalized() ;

		m_normals.col(i) = m_normalNoise.col(i).cast<float>() ;
		m_positions.col(i) = (particles.centers().col(pid) + m_offsets.col(i) ).cast< float >() ;
	}

}

void Sampler::reassign() 
{
	
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

#pragma omp parallel 
	{
		BallSampler sampler ;

		#pragma omp for
		for(Index i = 0 ; i < n ; ++i) {
			
			Vec v ;
			Mat frame ; 
			tensor_view( particles.frames().col(i) ).get( frame ) ;
			
			Eigen::SelfAdjointEigenSolver<Mat> es( frame );
			const Vec ev = es.eigenvalues().array().max(0).sqrt() ;

			for( unsigned j = 0 ; j < nSamples ; ++j )
			{
				const unsigned idx = i*nSamples +  j ;
				m_particleIds[idx] = i ;
				
				sampler(v) ;
				m_offsets.col( idx ) = ( es.eigenvectors() * ev.asDiagonal() ) * v * 1.7 ;
				sampler.onSphere( v ) ;

				m_normalNoise.col( idx ) = v ;

				m_normals.col( idx ) = m_normalNoise.col( idx ).cast<float>() ; //FIXME
				m_positions.col(idx) = (particles.centers().col(i) + m_offsets.col(idx) ).cast< float >() ;
			}
		}
	}
}


} // d6
