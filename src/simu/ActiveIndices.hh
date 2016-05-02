#ifndef D6_ACTIVE_INDICES_HH
#define D6_ACTIVE_INDICES_HH

#include "utils/alg.hh"
#include "geo/MeshBase.hh"

namespace d6 {

struct Active {

	static const Index s_Inactive  ;

	Index nNodes ; //!< Number of active nodes
	Index offset ; //!< Offset in global simulation nodes array
	typename MeshType::Cells cells ;  //!< List of active cells

	std::vector< Index > indices    ; //!< mesh_index -> simu_index
	std::vector< Index > revIndices ; //!< simu_index -> mesh_index

	Active()	 : nNodes( 0 ), offset(0) {}

	void reset( Index totNodes )
	{
		offset = 0 ;
		nNodes = 0 ;
		cells.clear();
		revIndices.clear() ;
		indices.assign( totNodes, s_Inactive );
	}

	void computeRevIndices() ;
	void setOffset( const Index o ) ;

	Index count() const { return nNodes ; }
	Index origSize() const { return indices.size() ; }

	//! Conversion from field on whole mesh to values at active nodes
	template < typename Derived, typename DestDerived >
	void field2var( const FieldBase<Derived> &field, Eigen::DenseBase<DestDerived>& var, bool resize = true ) const
	{
		constexpr Index D = FieldBase<Derived>::D ;

		if( resize )
			var.derived().resize( (offset + D) * count() );

#pragma omp parallel for
		for( Index i = 0 ; i < nNodes ; ++ i) {
			const Index idx = revIndices[ i ] ;
			Segmenter<D, DestDerived>::segment( var.derived(), offset+i ) = field[ idx ] ;
		}
	}

	//! Conversion from values at active nodes to field on whole mesh (zero at missing nodes)
	template < typename Derived >
	void var2field( const DynVec& var,  FieldBase<Derived> &field ) const
	{
		constexpr Index D = FieldBase<Derived>::D ;

		field.set_zero();

#pragma omp parallel for
		for( Index i = 0 ; i < nNodes ; ++ i) {
			const Index idx = revIndices[ i ] ;
			field[ idx ] = Segmenter<D>::segment( var, offset+i ) ;
		}
	}



};


} //d6


#endif
