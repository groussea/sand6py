#ifndef D6_TENSOR_FIELD_HH
#define D6_TENSOR_FIELD_HH

#include "FieldBase.hh"
#include "FieldFuncs.hh"

namespace d6 {

template < typename MeshT >
struct FieldTraits< AbstractTensorField< MeshT > >  {

	typedef MeshT  MeshType ;
	typedef VecS ValueType ;
	static constexpr Index Dimension = SD ;
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

	template <typename Func>
	AbstractTensorField( const FieldFuncBase< Func, Base::D, MeshT > & func )
		: Base( func.mesh() )
	{
		Base::operator=( func );
	}
	template <typename Func>
	AbstractTensorField& operator=( const FieldFuncBase< Func, Base::D, MeshT > & func )
	{
		return Base::operator=( func );
	}

	DeviatoricPart< MeshT > deviatoricPart() const {
		return DeviatoricPart< MeshT >( *this ) ;
	}
	FieldTrace< MeshT > trace() const {
		return FieldTrace< MeshT >( *this ) ;
	}
	FieldNorm< d6::AbstractTensorField, MeshT > norm() const {
		return FieldNorm< d6::AbstractTensorField, MeshT >( *this ) ;
	}

	void add_sym_tensor( const Vec& x, Mat& tensor ) const ;
	void get_sym_tensor( const Vec& x, Mat& tensor ) const ;

};


} //d6

#endif

