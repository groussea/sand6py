#ifndef D6_TENSOR_FIELD_HH
#define D6_TENSOR_FIELD_HH

#include "FieldBase.hh"
#include "FieldFuncs.hh"

namespace d6 {

// Symmetric tensor

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

	D6_MAKE_FIELD_CONSTRUCTORS_AND_ASSIGNMENT_OPERATORS( AbstractTensorField )

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

// Skew symmetric tensor

template < typename ShapeFuncT >
struct FieldTraits< AbstractSkewTsField< ShapeFuncT > >  {

	typedef ShapeFuncT  ShapeFuncImpl ;
	typedef VecR ValueType ;
	static constexpr Index Dimension = RD ;
};

template < typename ShapeFuncT >
class AbstractSkewTsField : public FieldBase< AbstractSkewTsField< ShapeFuncT > >
{
	typedef FieldBase< AbstractSkewTsField > Base ;

public:
	typedef ShapeFuncBase< ShapeFuncT > ShapeFuncType ;
	typedef typename ShapeFuncType::Location Location ;

	D6_MAKE_FIELD_CONSTRUCTORS_AND_ASSIGNMENT_OPERATORS( AbstractSkewTsField )

	void get_spi_tensor( const Location& x, Mat& tensor ) const ;
	void add_spi_tensor( const Location& x, Mat& tensor ) const ;

};


} //d6

#endif

