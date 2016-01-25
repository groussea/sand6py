#ifndef D6_VECTOR_FIELD_HH
#define D6_VECTOR_FIELD_HH

#include "FieldBase.hh"
#include "FieldFuncs.hh"

#include "Tensor.hh"

namespace d6 {

template < typename MeshT >
struct FieldTraits< AbstractVectorField< MeshT > >  {

	typedef MeshT  MeshType ;
	typedef Vec	   ValueType ;
	static constexpr Index Dimension = WD ;
};

template < typename MeshT >
class AbstractVectorField : public FieldBase< AbstractVectorField< MeshT > >
{
	typedef FieldTraits< AbstractVectorField > Traits ;
	typedef MeshBase< typename Traits::MeshType > MeshType ;

	typedef FieldBase< AbstractVectorField > Base ;

public:
	explicit AbstractVectorField( const MeshType& mesh )
		: Base( mesh )
	{

	}
	template <typename Func>
	AbstractVectorField( const FieldFuncBase< Func, Base::D, MeshT > & func )
		: Base( func.mesh() )
	{
		Base::operator=( func );
	}
	template <typename Func>
	AbstractVectorField& operator=( const FieldFuncBase< Func, Base::D, MeshT > & func )
	{
		return Base::operator=( func );
	}

	FieldNorm< d6::AbstractVectorField, MeshT > norm() const {
		return FieldNorm< d6::AbstractVectorField, MeshT >( *this ) ;
	}

#if (D6_DIM==3)
	void get_spi_tensor( const Vec& x, Mat& tensor ) const ;
	void add_spi_tensor( const Vec& x, Mat& tensor ) const ;
#endif

};


} //d6

#endif

