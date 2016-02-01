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

	explicit AbstractVectorField( const ShapeFuncType& shape )
		: Base( shape )
	{}
	explicit AbstractVectorField( const typename ShapeFuncType::DOFDefinition& mesh )
		: Base( ShapeFuncT( mesh ) )
	{}

	template <typename Func>
	AbstractVectorField( const FieldFuncBase< Func, Base::D, ShapeFuncT > & func )
		: Base( func.shape() )
	{
		Base::operator=( func );
	}
	template <typename Func>
	AbstractVectorField& operator=( const FieldFuncBase< Func, Base::D, ShapeFuncT > & func )
	{
		return Base::operator=( func );
	}

	FieldNorm< d6::AbstractVectorField, ShapeFuncT > norm() const {
		return FieldNorm< d6::AbstractVectorField, ShapeFuncT >( *this ) ;
	}

#if (D6_DIM==3)
	void get_spi_tensor( const Location& x, Mat& tensor ) const ;
	void add_spi_tensor( const Location& x, Mat& tensor ) const ;
#endif

};


} //d6

#endif

