#ifndef D6_TENSOR_FIELD_HH
#define D6_TENSOR_FIELD_HH

#include "FieldBase.hh"

namespace d6 {

template < typename MeshT >
class TensorField ;

template < typename MeshT >
struct FieldTraits< TensorField< MeshT > >  {

	typedef MeshT  MeshType ;
	typedef Vec6 ValueType ;
	static constexpr Index Dimension = 6 ;
};

template < typename MeshT >
class TensorField : public FieldBase< TensorField< MeshT > >
{
	typedef FieldTraits< TensorField > Traits ;
	typedef MeshBase< typename Traits::MeshType > MeshType ;

	typedef FieldBase< TensorField > Base ;

public:
	TensorField( const MeshType& mesh )
		: Base( mesh )
	{

	}

};


} //d6

#endif

