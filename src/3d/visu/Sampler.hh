#ifndef D6_SAMPLER_HH
#define D6_SAMPLER_HH

#include "utils/alg.hh"

#include <vector>

namespace d6 {

class Offline ;

class Sampler {

public:

	enum Mode {
		Normal,
		VelocityCut,
		Discs
	};

	explicit Sampler( const Offline& offline )
	: m_offline( offline ), m_particlesCount( 0 ), m_mode( Normal )
	{}



	void sampleParticles(unsigned nSamples ) ;

	void reassign() ;
	void move( ) ;

	size_t count() const { return m_particleIds.size() ; }

	const std::vector< unsigned >& particleIds() const { return m_particleIds ; }
	const DynMat3& offsets() const { return m_offsets ; }
	const Eigen::Matrix3Xf& noise() const { return m_normalNoise ; }

	const Eigen::Matrix3Xf& normals() const { return m_normals ; }
	const Eigen::Matrix3Xf& positions() const { return m_positions ; }
	const Eigen::VectorXf& visibility() const { return m_visibility ; }

	void setMode( const Mode mode ) {
		m_mode = mode ;
	}
	Mode mode() const {
		return m_mode ;
	}

	void compute_absolute( ) ;

private:

	const Offline& m_offline ;
	unsigned m_particlesCount ;

	Mode m_mode ;

	std::vector< unsigned > m_particleIds ;
	DynMat3 m_offsets ;
	Eigen::Matrix3Xf m_normalNoise ;

	Eigen::Matrix3Xf m_positions ;
	Eigen::Matrix3Xf m_normals ;
	Eigen::VectorXf m_visibility ;

	DynMat3 m_predPos ;

} ;

} // d6

#endif
