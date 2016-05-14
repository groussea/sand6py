#ifndef D6_FIELD_BASE_IMPL_HH
#define D6_FIELD_BASE_IMPL_HH

#include "FieldBase.hh"

#include "ShapeFunctionBase.hh"
#include "ScalarField.hh"

namespace d6 {


template< typename Derived >
void FieldBase< Derived >::add_at( const Location& x, const ValueType& val )
{
	typename ShapeFuncType::Interpolation itp ;
	m_shape.interpolate( x, itp );
	add_at( itp, val ) ;
}

template< typename Derived >
void FieldBase< Derived >::add_at( const typename ShapeFuncType::Interpolation &itp, const ValueType& val )
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
template< typename Derived >
Derived& FieldBase< Derived >::divide_by_positive(const ScalarField &field, Scalar min) {
	div_compwise< D >( m_data, field.flatten().cwiseMax(min) ) ;
	return derived() ;
}

template< typename Derived >
template< typename WorldFunc >
void FieldBase< Derived >::set_free( const WorldFunc&f )
{
	typedef typename ShapeFuncImpl::MeshType MeshType ;
	const MeshType& mesh = shape().mesh() ;

	Location dst_loc ;
	typename ShapeFuncType::NodeList cellNodes ;

	for( typename MeshType::CellIterator it = mesh.cellBegin() ;
		 it != mesh.cellEnd() ; ++it )
	{
		dst_loc.cell = *it ;
		shape().list_nodes( dst_loc, cellNodes );

		for( Index k = 0 ; k < ShapeFuncType::NI ; ++k ) {
			shape().locate_dof( dst_loc, k ) ;
			segment( cellNodes[k] ) = f( mesh.pos( dst_loc ) )  ;

		}
	}
}

// TODO fix reduduncy with/ non-free

template< typename Derived >
template< typename WorldFunc >
void FieldBase< Derived >::integrate_free( const WorldFunc&f )
{
	Location dst_loc ;

	for( auto qpIt = m_shape.qpBegin() ; qpIt != m_shape.qpEnd() ; ++qpIt ) {
		qpIt.locate( dst_loc ) ;
		add_at( dst_loc, qpIt.weight() * f( qpIt.pos() ) ) ;
	}
}

template< typename Derived >
template< typename WorldFunc >
void FieldBase< Derived >::interpolate_free( const WorldFunc&f )
{
	ScalarField volumes( shape() );
	shape().compute_lumped_mass( volumes.flatten() );

	set_zero() ;

	integrate_free( f ) ;

	// TODO : solve with consistent mass matrix ?
	divide_by_positive( volumes ) ;
}

template< typename Derived >
template< typename Func, typename OtherShape >
typename std::enable_if< OtherShape::is_mesh_based>::type
FieldBase< Derived >::integrate( const FieldFuncBase< Func, D, OtherShape > &f )
{
	typename OtherShape::Location src_loc ;
	Location dst_loc ;

	for( auto qpIt = m_shape.qpBegin() ; qpIt != m_shape.qpEnd() ; ++qpIt ) {
		qpIt.locate( dst_loc ) ;
		f.shape().derived().locate( qpIt.pos(), src_loc ) ;
		add_at( dst_loc, qpIt.weight() * f.eval_at( src_loc ) ) ;
	}
}

template< typename Derived >
template< typename Func, typename OtherShape >
typename std::enable_if<!OtherShape::is_mesh_based>::type
FieldBase< Derived >::integrate( const FieldFuncBase< Func, D, OtherShape > &f )
{
	typename OtherShape::Location src_loc ;
	Location dst_loc ;

	for( auto qpIt = f.shape().qpBegin() ; qpIt != f.shape().qpEnd() ; ++qpIt ) {
		qpIt.locate( src_loc ) ;
		m_shape.locate( qpIt.pos(), dst_loc ) ;
		add_at( dst_loc, qpIt.weight() * f.eval_at( src_loc ) ) ;
	}

}

template< typename Derived >
template< typename Func, typename OtherShape >
typename std::enable_if< FieldBase<Derived>::ShapeFuncType::is_mesh_based || OtherShape::is_mesh_based, Derived& >::type
FieldBase< Derived >::from_interpolation( const FieldFuncBase< Func, D, OtherShape > &f )
{
	ScalarField volumes( shape() );
	shape().compute_lumped_mass( volumes.flatten() );

	set_zero() ;

	integrate( f ) ;

	// TODO : solve with consistent mass matrix ?
	divide_by_positive( volumes ) ;

	return derived();
}

} //ns d6

#endif
