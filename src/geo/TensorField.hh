#ifndef D6_TENSOR_FIELD_HH
#define D6_TENSOR_FIELD_HH

#include "FieldBase.hh"

namespace d6 {

template < typename MeshT >
struct FieldTraits< AbstractTensorField< MeshT > >  {

	typedef MeshT  MeshType ;
	typedef Vec6 ValueType ;
	static constexpr Index Dimension = 6 ;
};

template < typename MeshT >
class AbstractTensorField : public FieldBase< AbstractTensorField< MeshT > >
{
	typedef FieldTraits< AbstractTensorField > Traits ;
	typedef MeshBase< typename Traits::MeshType > MeshType ;

	typedef FieldBase< AbstractTensorField > Base ;

public:
	explicit AbstractTensorField( const MeshType& mesh )
		: Base( mesh )
	{

	}

};


} //d6

#endif

