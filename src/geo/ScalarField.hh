#ifndef D6_SCALAR_FIELD_HH
#define D6_SCALAR_FIELD_HH

#include "FieldBase.hh"

namespace d6 {

template < typename MeshT >
class ScalarField ;

template < typename MeshT >
struct FieldTraits< ScalarField< MeshT > >  {

	typedef MeshT  MeshType ;
	typedef Scalar ValueType ;
	static constexpr Index Dimension = 1 ;
};

template < typename MeshT >
class ScalarField : public FieldBase< ScalarField< MeshT > >
{
	typedef FieldTraits< ScalarField > Traits ;
	typedef MeshBase< typename Traits::MeshType > MeshType ;

	typedef FieldBase< ScalarField > Base ;

public:
	ScalarField( const MeshType& mesh )
		: Base( mesh )
	{

	}

};


} //d6

#endif
