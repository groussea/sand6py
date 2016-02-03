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
	const typename VectorField::ShapeFuncImpl shape( m_mesh ) ;
	constexpr Index NV = VectorField::ShapeFuncType::NI ;

	Eigen::Matrix<float, WD, Eigen::Dynamic> vertices( WD, shape.nDOF() ) ;
	vertices.setZero() ;

	Eigen::Matrix<int, NV+1, Eigen::Dynamic > nodeIndices( NV+1, m_mesh.nCells() ) ;
	nodeIndices.row(0).setConstant( NV ) ;

	const Index elemType =
		#if D6_DIM == 3
		MeshType::NV == 4 ? /*VTK_TETRA*/ 10 : /*VTK_VOXEL*/ 11
		#else
		MeshType::NV == 3 ? /*VTK_TRIANGLE*/ 5 : /*VTK_PIXEL*/ 8
		#endif
			;

	const Eigen::VectorXi cellTypes =
		Eigen::VectorXi::Constant( m_mesh.nCells(), elemType ) ;


	typename MeshType::CellGeo cellGeo ;
	typename  VectorField::ShapeFuncType::NodeList cellNodes ;
	typename  VectorField::ShapeFuncType::Location loc ;

	for( typename MeshType::CellIterator it = m_mesh.cellBegin() ; it != m_mesh.cellEnd() ; ++it )
	{
		m_mesh.get_geo( *it, cellGeo ) ;
		loc.cell = *it ;
		shape.list_nodes( loc, cellNodes );

		nodeIndices.block< NV, 1 >( 1, it.index() ) = cellNodes ;

		for( int k = 0 ; k < NV ; ++k ) {
			vertices.col( cellNodes[k] ).head<WD>() = cellGeo.vertex( k ).cast< float >() ;
		}
	}

	vtk << "DATASET UNSTRUCTURED_GRID\n" ;
	vtk << "POINTS " << vertices.cols() << " float\n" ;
	write( vtk, vertices.data(), WD, vertices.cols() ) ;
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
