#ifndef D6_SAMPLER_HH
#define D6_SAMPLER_HH

#include "utils/alg.hh"

#include <vector>

namespace d6 {

class Offline ;

class Sampler {

public:
	explicit Sampler( const Offline& offline )
	: m_offline( offline )
	{}

	void sampleParticles( unsigned nSamples ) ;

	void reassign() ;
	void updateOffsets( const Scalar dt ) ;

	size_t count() const { return m_particleIds.size() ; }

	const std::vector< unsigned >& particleIds() const { return m_particleIds ; }
	const DynMat3& offsets() const { return m_offsets ; }

	const Eigen::Matrix3Xf& normals() const { return m_normals ; } 
	const Eigen::Matrix3Xf& positions() const { return m_positions ; } 
	

private:

	const Offline& m_offline ;

	std::vector< unsigned > m_particleIds ;
	DynMat3 m_offsets ;
	DynMat3 m_normalNoise ;
	
	Eigen::Matrix3Xf m_positions ;
	Eigen::Matrix3Xf m_normals ;

} ;

} // d6

#endif

