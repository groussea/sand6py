#ifndef D6_FIELD_FUNCS_IMPL_HPP
#define D6_FIELD_FUNCS_IMPL_HPP

#include "FieldFuncBase.hh"

namespace d6 {

template< typename Derived >
struct FuncTraits {
};

template < typename Derived, Index D,
		   typename MeshT = typename FieldTraits< typename FuncTraits<Derived>::FieldType >::MeshType >
struct UnaryFieldFunc : public FieldFuncBase< Derived, D, MeshT >
{
	typedef FieldFuncBase< Derived, D, MeshT > Base ;
	typedef FuncTraits< Derived > Traits;
	typedef typename Traits::FieldType FieldType;

	explicit UnaryFieldFunc( const FieldType & field )
		: Base( field.mesh() ), m_field(field)
	{}

	Index size( ) const {
		return m_field.size() ;
	}

protected:
	const FieldType& m_field  ;
};

// Trace

template <typename MeshT> struct FieldTrace ;

template< typename MeshT >
struct FuncTraits< FieldTrace< MeshT > > {
	typedef AbstractTensorField< MeshT > FieldType ;
};

template <typename MeshT>
struct FieldTrace : public UnaryFieldFunc<FieldTrace<MeshT>, 1>
{
	typedef UnaryFieldFunc<FieldTrace<MeshT>, 1> Base ;
	typedef typename Base::FieldType FieldType ;

	using Base::m_field ;
	explicit FieldTrace( const FieldType & field ) : Base(field) {}

	void eval_at_node( Index i, typename Base::Seg v ) const
	{
		v = m_field[i][0] ;
	}

} ;

// Deviatoric part

template <typename MeshT> struct DeviatoricPart ;

template< typename MeshT >
struct FuncTraits< DeviatoricPart< MeshT > > {
	typedef AbstractTensorField< MeshT > FieldType ;
};

template <typename MeshT>
struct DeviatoricPart : public UnaryFieldFunc<DeviatoricPart<MeshT>, SD>
{
	typedef UnaryFieldFunc<DeviatoricPart<MeshT>, SD> Base ;
	typedef typename Base::FieldType FieldType ;

	using Base::m_field ;
	explicit DeviatoricPart( const FieldType & field ) : Base(field) {}

	void eval_at_node( Index i, typename Base::Seg v ) const
	{
		v = m_field[i] ;
		v[0] = 0 ;
	}

} ;

// Norm

template <template <typename> class Field, typename MeshT >
struct FieldNorm ;

template <template <typename> class Field, typename MeshT >
struct FuncTraits< FieldNorm< Field, MeshT > > {
	typedef Field< MeshT > FieldType ;
};

template <template <typename> class Field, typename MeshT >
struct FieldNorm : public UnaryFieldFunc< FieldNorm<Field, MeshT>, 1>
{
	typedef UnaryFieldFunc<FieldNorm<Field, MeshT>, 1> Base ;
	typedef typename Base::FieldType FieldType ;

	using Base::m_field ;
	explicit FieldNorm( const FieldType & field ) : Base(field) {}

	void eval_at_node( Index i, typename Base::Seg v ) const
	{
		v = m_field[i].norm() ;
	}

} ;


} //d6

#endif
