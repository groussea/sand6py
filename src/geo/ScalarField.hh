#ifndef D6_SCALAR_FIELD_HH
#define D6_SCALAR_FIELD_HH

#include "FieldBase.hh"

namespace d6 {

template < typename MeshT >
struct FieldTraits< AbstractScalarField< MeshT > >  {

	typedef MeshT  MeshType ;
	typedef Scalar ValueType ;
	static constexpr Index Dimension = 1 ;
};

template < typename MeshT >
class AbstractScalarField : public FieldBase< AbstractScalarField< MeshT > >
{
	typedef FieldTraits< AbstractScalarField > Traits ;
	typedef MeshBase< typename Traits::MeshType > MeshType ;

	typedef FieldBase< AbstractScalarField > Base ;

public:
	AbstractScalarField( const MeshType& mesh )
		: Base( mesh )
	{

	}

};


} //d6

#endif
