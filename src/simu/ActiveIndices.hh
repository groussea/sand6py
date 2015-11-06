#ifndef D6_ACTIVE_INDICES_HH
#define D6_ACTIVE_INDICES_HH

#include "utils/alg.hh"
#include "geo/MeshBase.hh"

namespace d6 {

struct Active {

	static const Index s_Inactive  ;

	Index nNodes ;
	Index offset ;
	typename MeshType::Cells cells ;

	std::vector< Index > indices    ;
	std::vector< Index > revIndices ;

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

	template < typename Derived >
	void field2var( const FieldBase<Derived> &field, DynVec & var, bool resize = true ) const ;
	template < typename Derived >
	void var2field( const DynVec & var, FieldBase<Derived> &field ) const ;

};


} //d6


#endif
