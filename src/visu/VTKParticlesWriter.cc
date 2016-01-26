#include "VTKParticlesWriter.hh"

#include "geo/Particles.hh"

#include "utils/File.hh"
#include "utils/Log.hh"

namespace d6 {

VTKParticlesWriter::VTKParticlesWriter(const char *base_dir, const Particles &particles)
	: VTKWriter( base_dir ), m_particles(particles)
{

}

void VTKParticlesWriter::writeMesh( File &vtk ) const
{
	vtk << "DATASET POLYDATA\n" ;
	vtk << "POINTS " << m_particles.count() << " float\n" ;
	write( vtk, m_particles.centers().data(), WD, m_particles.count() ) ;
}

size_t VTKParticlesWriter::nDataPoints() const {
	return m_particles.count() ;
}

template< typename Derived >
bool VTKParticlesWriter::dump( const char *name, const Eigen::MatrixBase< Derived > &data )
{
	if( !m_file.is_open() ) {
		Log::Error() << " VTKParticlesWriter: should call startFile() before dumping data " << std::endl ;
		return false ;
	}

	writeAttribute( name, data.derived().data(), Derived::RowsAtCompileTime ) ;

	return true ;
}

bool VTKParticlesWriter::dump_all( )
{
	return
			dump( Volumes ) &&
			dump( Velocities ) &&
			dump( Frames ) &&
			dump( Orientations ) ;

	return true ;
}

bool VTKParticlesWriter::dump( Quantity quantity) {
	switch( quantity ) {
	case Volumes:
		return dump( "volumes", m_particles.volumes() ) ;
	case Velocities:
		return dump( "velocities", m_particles.velocities() ) ;
	case Frames:
		return dump( "frames", m_particles.frames() ) ;
	case Orientations:
		return dump( "orient", m_particles.orient() ) ;
	}
	return false ;
}

}
