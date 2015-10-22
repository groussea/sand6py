#ifndef D6_VECTOR_FIELD_HH
#define D6_VECTOR_FIELD_HH

#include "FieldBase.hh"

namespace d6 {

template < typename MeshT >
struct FieldTraits< AbstractVectorField< MeshT > >  {

	typedef MeshT  MeshType ;
	typedef Vec	   ValueType ;
	static constexpr Index Dimension = 3 ;
};

template < typename MeshT >
class AbstractVectorField : public FieldBase< AbstractVectorField< MeshT > >
{
	typedef FieldTraits< AbstractVectorField > Traits ;
	typedef MeshBase< typename Traits::MeshType > MeshType ;

	typedef FieldBase< AbstractVectorField > Base ;

public:
	explicit AbstractVectorField( const MeshType& mesh )
		: Base( mesh )
	{

	}

};


} //d6

#endif

