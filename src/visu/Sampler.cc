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

void Sampler::move( const Scalar dt )
{
	// offset += DU*offset * st

	const Particles& particles = m_offline.particles() ;
	const Phase& grains = m_offline.grains() ;

	const Index n = count() ;
	
#pragma omp parallel for
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
	}

	compute_absolute() ;
}

void Sampler::compute_absolute( )
{
	// offset += DU*offset * st

	const Scalar noise = 1. ;

	const Particles& particles = m_offline.particles() ;
	const Phase& grains = m_offline.grains() ;

	const Index n = count() ;
	
#pragma omp parallel for
	for(Index i = 0 ; i < n ; ++i) {
		const unsigned pid = m_particleIds[i] ;
		const Vec p0 = particles.centers().col( pid ) ;

		const Vec pos = (p0 + m_offsets.col(i) ) ;
		m_positions.col(i) = pos.cast< float >() ;
		
		Eigen::Vector3f grad_phi = grains.grad_phi( pos ).cast<float>() ;
		grad_phi += 1.e-2 * Eigen::Vector3f( 0, 0, -1 ) ;
		const float gn = grad_phi.norm() ;
		if( gn > 1.e-6 )
			grad_phi /= gn ;
		else
			grad_phi = Eigen::Vector3f( 0, 0, -1 ) ;


		Eigen::Quaternionf rot( Eigen::AngleAxisf( noise, m_normalNoise.col(i).cast< float >() ) ) ;
		m_normals.col(i) =  - ( rot * grad_phi ).normalized() ;
		
		m_visibility(i) = std::max( 0., std::min( 1., 1.5 - 1.5*grains.fraction(pos) ) ) ;
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
	m_visibility.resize( n * nSamples ) ;

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

			}
		}
	}

	compute_absolute() ;
}


} // d6
