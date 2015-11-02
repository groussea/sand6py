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
	explicit AbstractScalarField( const MeshType& mesh )
		: Base( mesh )
	{

	}

	template <typename Func>
	AbstractScalarField( const FieldFuncBase< Func, Base::D, MeshT > & func )
		: Base( func.mesh() )
	{
		Base::operator=( func );
	}
	template <typename Func>
	AbstractScalarField& operator= ( const FieldFuncBase< Func, Base::D, MeshT > & func )
	{
		return Base::operator=( func );
	}

};


} //d6

#endif
