#include "VTKFieldWriter.hh"

#include "geo/ScalarField.hh"
#include "geo/VectorField.hh"
#include "geo/TensorField.hh"

#include "geo/MeshImpl.hh"

#include "utils/File.hh"
#include "utils/Log.hh"

namespace d6 {

VTKFieldWriter::VTKFieldWriter( const char* base_dir, const MeshType& mesh )
	: VTKWriter( base_dir ), m_mesh( mesh )
{
}

void VTKFieldWriter::writeMesh( File &vtk ) const
{

	Eigen::Matrix3Xf vertices( 3, m_mesh.nNodes() ) ;

	Eigen::Matrix<int, MeshType::NV+1, Eigen::Dynamic > nodeIndices( MeshType::NV+1, m_mesh.nCells() ) ;
	nodeIndices.row(0).setConstant( MeshType::NV ) ;

	const Eigen::VectorXi cellTypes = Eigen::VectorXi::Constant( m_mesh.nCells(), /*VTK_VOXEL*/ 11 ) ;

	typename MeshType::CellGeo cellGeo ;
	typename MeshType::NodeList cellNodes ;

	for( typename MeshType::CellIterator it = m_mesh.cellBegin() ; it != m_mesh.cellEnd() ; ++it )
	{
		m_mesh.get_geo( *it, cellGeo ) ;
		m_mesh.list_nodes( *it, cellNodes );

		nodeIndices.block< MeshType::NV, 1 >( 1, it.index() ) = cellNodes ;

		for( int k = 0 ; k < MeshType::NV ; ++k ) {
			vertices.col( cellNodes[k] ) = cellGeo.vertex( k ).cast< float >() ;
		}
	}

	vtk << "DATASET UNSTRUCTURED_GRID\n" ;
	vtk << "POINTS " << vertices.cols() << " float\n" ;
	write( vtk, vertices.data(), 3, vertices.cols() ) ;
	vtk << "CELLS " << nodeIndices.cols() << " " << nodeIndices.size() << "\n";
	write( vtk, nodeIndices.data(), 1, nodeIndices.size() ) ;
	vtk << "CELL_TYPES " << nodeIndices.cols() << "\n";
	write( vtk, cellTypes.data(), 1, nodeIndices.cols() ) ;
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
