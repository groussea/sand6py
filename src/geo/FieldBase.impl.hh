#ifndef D6_FIELD_BASE_IMPL_HH
#define D6_FIELD_BASE_IMPL_HH

#include "FieldBase.hh"

namespace d6 {

template< typename Derived >
void FieldBase< Derived >::eval_at( const Vec& x, ValueType& res ) const
{
	typename MeshType::Location loc ;
	m_mesh.locate( x, loc );

	d6::set_zero( res ) ;
	for( Index k = 0 ; k < loc.nodes.rows() ; ++k ) {
		res += loc.coeffs[k] * segment( loc.nodes[k] ) ;
	}

}

template< typename Derived >
void FieldBase< Derived >::add_at( const Vec& x, const ValueType& val )
{
	typename MeshType::Location loc ;
	m_mesh.locate( x, loc );

	for( Index k = 0 ; k < loc.nodes.rows() ; ++k ) {
		segment( loc.nodes[k] ) += loc.coeffs[k] * val ;
	}

}

template< typename Derived >
void FieldBase< Derived >::set_zero() {
	m_data.setZero() ;
}

template< typename Derived >
void FieldBase< Derived >::set_constant(const ValueType &val) {
#pragma omp parallel for
	for( Index i = 0 ; i < m_size ; ++i ) {
		segment(i) = val ;
	}
}

} //ns d6

#endif
