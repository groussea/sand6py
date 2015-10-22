#include "VTKFieldWriter.hh"

#include "geo/ScalarField.hh"
#include "geo/VectorField.hh"
#include "geo/TensorField.hh"

#include "geo/Grid.hh"

#include "utils/File.hh"
#include "utils/Log.hh"

namespace d6 {

VTKFieldWriter::VTKFieldWriter( const char* base_dir, const MeshType& mesh )
	: VTKWriter( base_dir ), m_mesh( mesh )
{
}

void VTKFieldWriter::writeMesh( File &vtk ) const
{
	vtk << "DATASET STRUCTURED_POINTS\n" ;
	const Grid& g = m_mesh.derived() ;
	vtk << "DIMENSIONS " << (g.dim() + Vec3i::Ones()).transpose() << "\n" ;
	vtk << "ORIGIN " << Vec::Zero().transpose() << "\n" ;
	vtk << "SPACING " << g.dx().transpose() << "\n" ;
}

size_t VTKFieldWriter::nDataPoints() const
{
	return m_mesh.nNodes() ;
}

template< typename Derived >
bool VTKFieldWriter::dump( const char* name, const FieldBase< Derived >& field )
{
	if( !m_file.is_open() ) {
		Log::Error() << " VTKParticlesWriter: should call startFile() before dumping data " << std::endl ;
		return false ;
	}

	writeAttribute( name, field.flatten().data(), FieldBase< Derived >::D ) ;

	return true ;
}

template bool VTKFieldWriter::dump( const char*, const FieldBase< ScalarField >& ) ;
template bool VTKFieldWriter::dump( const char*, const FieldBase< VectorField >& ) ;
template bool VTKFieldWriter::dump( const char*, const FieldBase< TensorField >& ) ;

} //d6
