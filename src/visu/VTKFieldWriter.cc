#include "VTKFieldWriter.hh"

#include "geo/ScalarField.hh"
#include "geo/VectorField.hh"
#include "geo/TensorField.hh"

#include "geo/MeshImpl.hh"

#include "utils/File.hh"
#include "utils/Log.hh"

namespace d6 {

template <typename ShapeFuncT>
VTKFieldWriter<ShapeFuncT>::VTKFieldWriter(const char* base_dir, const ShapeFuncT &shape )
	: VTKWriter( base_dir ), m_shape( shape )
{
}

template <typename ShapeFuncT>
void VTKFieldWriter<ShapeFuncT>::writeMesh( File &vtk ) const
{
	constexpr Index NV = ShapeFuncT::NI ;
	typedef typename ShapeFuncT::MeshType MeshType ;
	const MeshType& mesh = m_shape.mesh() ;

	Eigen::Matrix<float, WD, Eigen::Dynamic> vertices( WD, m_shape.nDOF() ) ;
	vertices.setZero() ;

	Eigen::Matrix<int, NV+1, Eigen::Dynamic > nodeIndices( NV+1, mesh.nCells() ) ;
	nodeIndices.row(0).setConstant( NV ) ;

	const Index elemType =
		#if D6_DIM == 3
		MeshType::NV == 4 ? /*VTK_TETRA*/ 10 : /*VTK_VOXEL*/ 11
		#else
		MeshType::NV == 3 ? /*VTK_TRIANGLE*/ 5 : /*VTK_PIXEL*/ 8
		#endif
			;

	const Eigen::VectorXi cellTypes =
		Eigen::VectorXi::Constant( mesh.nCells(), elemType ) ;


	typename MeshType::CellGeo cellGeo ;
	typename ShapeFuncT::NodeList cellNodes ;
	typename ShapeFuncT::Location loc ;

	// FXIME if shape nodes do not coincide with mesh
	for( typename MeshType::CellIterator it = mesh.cellBegin() ; it != mesh.cellEnd() ; ++it )
	{
		mesh.get_geo( *it, cellGeo ) ;
		loc.cell = *it ;
		m_shape.list_nodes( loc, cellNodes );

		nodeIndices.template block< NV, 1 >( 1, it.index() ) = cellNodes ;

		for( int k = 0 ; k < NV ; ++k ) {
			vertices.col( cellNodes[k] ).template head<WD>() = cellGeo.vertex( k ).cast< float >() ;
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

template <typename ShapeFuncT>
size_t VTKFieldWriter<ShapeFuncT>::nDataPoints() const
{
	return m_shape.nDOF() ;
}

template <typename ShapeFuncT>
template< typename Derived >
bool VTKFieldWriter<ShapeFuncT>::dump( const char* name, const FieldBase< Derived >& field )
{
	if( !m_file.is_open() ) {
		Log::Error() << " VTKParticlesWriter: should call startFile() before dumping data " << std::endl ;
		return false ;
	}

	writeAttribute( name, field.flatten().data(), FieldBase< Derived >::D ) ;

	return true ;
}

template class VTKFieldWriter<Linear<MeshImpl>> ;
template  bool VTKFieldWriter<Linear<MeshImpl>>::dump( const char*, const FieldBase< AbstractScalarField<Linear<MeshImpl> > >& ) ;
template  bool VTKFieldWriter<Linear<MeshImpl>>::dump( const char*, const FieldBase< AbstractVectorField<Linear<MeshImpl> > >& ) ;
template  bool VTKFieldWriter<Linear<MeshImpl>>::dump( const char*, const FieldBase< AbstractTensorField<Linear<MeshImpl> > >& ) ;

template class VTKFieldWriter<DGLinear<MeshImpl>> ;
template  bool VTKFieldWriter<DGLinear<MeshImpl>>::dump( const char*, const FieldBase< AbstractScalarField<DGLinear<MeshImpl> > >& ) ;
template  bool VTKFieldWriter<DGLinear<MeshImpl>>::dump( const char*, const FieldBase< AbstractTensorField<DGLinear<MeshImpl> > >& ) ;

} //d6
