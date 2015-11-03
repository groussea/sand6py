#ifndef D6_FIELD_BASE_IMPL_HH
#define D6_FIELD_BASE_IMPL_HH

#include "FieldBase.hh"

#include "MeshBase.hh"
#include "ScalarField.hh"

namespace d6 {

template< typename Derived >
void FieldBase< Derived >::eval_at( const Vec& x, ValueType& res ) const
{
	typename MeshType::Interpolation itp ;
	m_mesh.interpolate( x, itp );

	d6::set_zero( res ) ;
	for( Index k = 0 ; k < itp.nodes.rows() ; ++k ) {
		res += itp.coeffs[k] * segment( itp.nodes[k] ) ;
	}

}

template< typename Derived >
void FieldBase< Derived >::add_at( const Vec& x, const ValueType& val )
{
	typename MeshType::Interpolation itp ;
	m_mesh.interpolate( x, itp );
	add_at( itp, val ) ;
}

template< typename Derived >
void FieldBase< Derived >::add_at( const typename MeshType::Interpolation &itp, const ValueType& val )
{
	for( Index k = 0 ; k < itp.nodes.rows() ; ++k ) {
		segment( itp.nodes[k] ) += itp.coeffs[k] * val ;
	}

}

template< typename Derived >
Derived& FieldBase< Derived >::set_zero() {
	m_data.setZero() ;
	return derived() ;
}

template< typename Derived >
Derived& FieldBase< Derived >::set_constant(const ValueType &val) {
#pragma omp parallel for
	for( Index i = 0 ; i < m_size ; ++i ) {
		segment(i) = val ;
	}
	return derived() ;
}
template< typename Derived >
Derived& FieldBase< Derived >::multiply_by(const ScalarField &field) {
	mul_compwise< D >( m_data, field.flatten() ) ;
	return derived() ;
}
template< typename Derived >
Derived& FieldBase< Derived >::divide_by(const ScalarField &field) {
	div_compwise< D >( m_data, field.flatten() ) ;
	return derived() ;
}

} //ns d6

#endif
