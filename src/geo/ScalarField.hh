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

	explicit AbstractScalarField( const ShapeFuncType& shape )
		: Base( shape )
	{

	}
	explicit AbstractScalarField( const typename ShapeFuncType::DOFDefinition& mesh )
		: Base( ShapeFuncT( mesh ) )
	{}

	template <typename Func>
	AbstractScalarField( const FieldFuncBase< Func, Base::D, ShapeFuncT > & func )
		: Base( func.shape() )
	{
		Base::operator=( func );
	}
	template <typename Func>
	AbstractScalarField& operator= ( const FieldFuncBase< Func, Base::D, ShapeFuncT > & func )
	{
		return Base::operator=( func );
	}

	Vec grad_at( const typename ShapeFuncType::Location& loc ) const ;

#if (D6_DIM==2)
	void get_spi_tensor( const Location& x, Mat& tensor ) const ;
	void add_spi_tensor( const Location& x, Mat& tensor ) const ;
#endif

};


} //d6

#endif
