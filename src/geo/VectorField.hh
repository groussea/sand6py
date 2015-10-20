#ifndef D6_VECTOR_FIELD_HH
#define D6_VECTOR_FIELD_HH

#include "FieldBase.hh"

namespace d6 {

template < typename MeshT >
class VectorField ;

template < typename MeshT >
struct FieldTraits< VectorField< MeshT > >  {

	typedef MeshT  MeshType ;
	typedef Vec	   ValueType ;
	static constexpr Index Dimension = 3 ;
};

template < typename MeshT >
class VectorField : public FieldBase< VectorField< MeshT > >
{
	typedef FieldTraits< VectorField > Traits ;
	typedef MeshBase< typename Traits::MeshType > MeshType ;

	typedef FieldBase< VectorField > Base ;

public:
	VectorField( const MeshType& mesh )
		: Base( mesh )
	{

	}

};


} //d6

#endif

