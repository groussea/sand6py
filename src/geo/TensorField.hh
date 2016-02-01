#ifndef D6_TENSOR_FIELD_HH
#define D6_TENSOR_FIELD_HH

#include "FieldBase.hh"
#include "FieldFuncs.hh"

namespace d6 {

template < typename ShapeFuncT >
struct FieldTraits< AbstractTensorField< ShapeFuncT > >  {

	typedef ShapeFuncT  ShapeFuncImpl ;
	typedef VecS ValueType ;
	static constexpr Index Dimension = SD ;
};

template < typename ShapeFuncT >
class AbstractTensorField : public FieldBase< AbstractTensorField< ShapeFuncT > >
{
	typedef FieldBase< AbstractTensorField > Base ;

public:
	typedef ShapeFuncBase< ShapeFuncT > ShapeFuncType ;
	typedef typename ShapeFuncType::Location Location ;

	explicit AbstractTensorField( const ShapeFuncType& shape )
		: Base( shape )
	{
	}
	explicit AbstractTensorField( const typename ShapeFuncType::DOFDefinition& mesh )
		: Base( ShapeFuncT( mesh ) )
	{}

	template <typename Func>
	AbstractTensorField( const FieldFuncBase< Func, Base::D, ShapeFuncT > & func )
		: Base( func.shape() )
	{
		Base::operator=( func );
	}
	template <typename Func>
	AbstractTensorField& operator=( const FieldFuncBase< Func, Base::D, ShapeFuncT > & func )
	{
		return Base::operator=( func );
	}

	DeviatoricPart< ShapeFuncT > deviatoricPart() const {
		return DeviatoricPart< ShapeFuncT >( *this ) ;
	}
	FieldTrace< ShapeFuncT > trace() const {
		return FieldTrace< ShapeFuncT >( *this ) ;
	}
	FieldNorm< d6::AbstractTensorField, ShapeFuncT > norm() const {
		return FieldNorm< d6::AbstractTensorField, ShapeFuncT >( *this ) ;
	}

	void add_sym_tensor( const Location& x, Mat& tensor ) const ;
	void get_sym_tensor( const Location& x, Mat& tensor ) const ;

};


} //d6

#endif

