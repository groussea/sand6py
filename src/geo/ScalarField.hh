#ifndef D6_SCALAR_FIELD_HH
#define D6_SCALAR_FIELD_HH

#include "FieldBase.hh"

namespace d6 {

template < typename ShapeFuncT >
struct FieldTraits< AbstractScalarField< ShapeFuncT > >  {

	typedef ShapeFuncT ShapeFuncImpl ;
	typedef Scalar ValueType ;
	static constexpr Index Dimension = 1 ;
};

template < typename ShapeFuncT >
class AbstractScalarField : public FieldBase< AbstractScalarField< ShapeFuncT > >
{
	typedef FieldBase< AbstractScalarField > Base ;

public:
	typedef ShapeFuncBase< ShapeFuncT > ShapeFuncType ;
	typedef typename ShapeFuncType::Location Location ;

	D6_MAKE_FIELD_CONSTRUCTORS_AND_ASSIGNMENT_OPERATORS( AbstractScalarField )

	Vec grad_at( const typename ShapeFuncType::Location& loc ) const ;

};


} //d6

#endif
