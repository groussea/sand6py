#include "VTKParticlesWriter.hh"

#include "geo/Particles.hh"

#include "utils/File.hh"
#include "utils/Log.hh"

namespace d6 {

VTKParticlesWriter::VTKParticlesWriter(const char *base_dir, const Particles &particles)
	: VTKWriter( base_dir ), m_particles(particles)
{

}

void VTKParticlesWriter::writePoints( File &vtk ) const
{
	vtk << "DATASET POLYDATA\n" ;
	vtk << "POINTS " << m_particles.count() << " float\n" ;
	write( vtk, m_particles.centers().data(), 3, m_particles.count() ) ;
}

template< typename Derived >
bool VTKParticlesWriter::dump( unsigned frame, const char*name, const Eigen::MatrixBase< Derived > &data ) const
{
	File vtk ;
	if( !open( frame, name, vtk ) )
		return false ;

	writeHeader( vtk, name ) ;
	writePoints( vtk );

	vtk << "POINT_DATA " << m_particles.count() << "\n" ;
	writeDataHeader( vtk, Derived::RowsAtCompileTime, name ) ;
	write( vtk, data.derived().data(), Derived::RowsAtCompileTime, m_particles.count() ) ;

	return true ;
}

bool VTKParticlesWriter::dump_all( unsigned frame ) const
{
	const char *name = "particles" ;

	File vtk ;
	if( !open( frame, name, vtk ) )
		return false ;

	writeHeader( vtk, name ) ;
	writePoints( vtk );

	vtk << "POINT_DATA " << m_particles.count() << "\n" ;
	writeDataHeader( vtk, 1, "volume" ) ;
	write( vtk, m_particles.volumes().data(), 1, m_particles.count() ) ;
	writeDataHeader( vtk, 3, "velocity" ) ;
	write( vtk, m_particles.velocities().data(), 3, m_particles.count() ) ;
	writeDataHeader( vtk, 6, "frames" ) ;
	write( vtk, m_particles.frames().data(), 6, m_particles.count() ) ;
	writeDataHeader( vtk, 6, "orient" ) ;
	write( vtk, m_particles.orient().data(), 6, m_particles.count() ) ;

	return true ;
}

bool VTKParticlesWriter::dump(unsigned frame, Quantity quantity) const {
	switch( quantity ) {
	case Volumes:
		return dump( frame, "p_volumes", m_particles.volumes() ) ;
	case Velocities:
		return dump( frame, "p_velocities", m_particles.velocities() ) ;
	case Frames:
		return dump( frame, "p_frames", m_particles.frames() ) ;
	case Orientations:
		return dump( frame, "p_orient", m_particles.orient() ) ;
	}
	return false ;
}

}
