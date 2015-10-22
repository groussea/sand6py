#ifndef D6_VTK_PARTICLES_WRITER_HH
#define D6_VTK_PARTICLES_WRITER_HH

#include "VTKWriter.hh"
#include "utils/alg.hh"

namespace d6 {

class Particles ;

class VTKParticlesWriter : public VTKWriter
{

public:
	enum Quantity {
		Volumes,
		Velocities,
		Frames,
		Orientations
	} ;

	VTKParticlesWriter( const char* base_dir, const Particles& particles ) ;

	template< typename Derived >
	bool dump( unsigned frame, const char* name, const Eigen::MatrixBase< Derived > &data ) const ;

	bool dump( unsigned frame, Quantity quantity ) const ;

private:
	void writePoints(File &file) const ;

	const Particles& m_particles ;

};

} //d6

#endif
