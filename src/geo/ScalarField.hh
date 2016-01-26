#ifndef D6_SCALAR_FIELD_HH
#define D6_SCALAR_FIELD_HH

#include "FieldBase.hh"

namespace d6 {

template < typename MeshT >
struct FieldTraits< AbstractScalarField< MeshT > >  {

	typedef MeshT  MeshType ;
	typedef Scalar ValueType ;
	static constexpr Index Dimension = 1 ;
};

template < typename MeshT >
class AbstractScalarField : public FieldBase< AbstractScalarField< MeshT > >
{
	typedef FieldTraits< AbstractScalarField > Traits ;
	typedef MeshBase< typename Traits::MeshType > MeshType ;

	typedef FieldBase< AbstractScalarField > Base ;

public:
	explicit AbstractScalarField( const MeshType& mesh )
		: Base( mesh )
	{

	}

	// Naive local gradient
	Vec grad_at( const Vec& x) const ;

	template <typename Func>
	AbstractScalarField( const FieldFuncBase< Func, Base::D, MeshT > & func )
		: Base( func.mesh() )
	{
		Base::operator=( func );
	}
	template <typename Func>
	AbstractScalarField& operator= ( const FieldFuncBase< Func, Base::D, MeshT > & func )
	{
		return Base::operator=( func );
	}

#if (D6_DIM==2)
	void get_spi_tensor( const Vec& x, Mat& tensor ) const ;
	void add_spi_tensor( const Vec& x, Mat& tensor ) const ;
#endif

protected:
	using Base::m_mesh ;
};


} //d6

#endif
