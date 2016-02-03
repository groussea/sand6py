#ifndef D6_VECTOR_FIELD_HH
#define D6_VECTOR_FIELD_HH

#include "FieldBase.hh"
#include "FieldFuncs.hh"

#include "Tensor.hh"

namespace d6 {

template < typename ShapeFuncT >
struct FieldTraits< AbstractVectorField< ShapeFuncT > >  {

	typedef ShapeFuncT  ShapeFuncImpl ;
	typedef Vec	   ValueType ;
	static constexpr Index Dimension = WD ;
};

template < typename ShapeFuncT >
class AbstractVectorField : public FieldBase< AbstractVectorField< ShapeFuncT > >
{
	typedef FieldBase< AbstractVectorField > Base ;

public:
	typedef ShapeFuncBase< ShapeFuncT > ShapeFuncType ;
	typedef typename ShapeFuncType::Location Location ;

	D6_MAKE_FIELD_CONSTRUCTORS_AND_ASSIGNMENT_OPERATORS( AbstractVectorField )

	FieldNorm< d6::AbstractVectorField, ShapeFuncT > norm() const {
		return FieldNorm< d6::AbstractVectorField, ShapeFuncT >( *this ) ;
	}

};


} //d6

#endif

