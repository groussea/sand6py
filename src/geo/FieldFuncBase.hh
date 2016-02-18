#ifndef D6_FIELD_FUNC_HH
#define D6_FIELD_FUNC_HH

#include "utils/alg.hh"
#include "geo.fwd.hh"

namespace d6 {

template < typename Derived, Index D, typename ShapeFuncT, typename DestShapeFuncT >
struct Interpolation ;

template < typename Derived, Index D, typename ShapeFuncT >
struct FieldFuncBase
{
	typedef ShapeFuncBase< ShapeFuncT > ShapeFuncType ;
	typedef typename ShapeFuncType::Location Location ;

	typedef typename Segmenter< D >::ValueType ValueType ;

	explicit FieldFuncBase( const ShapeFuncType& shape ) : m_shape(shape.derived()) {}

	// Evaluation at arbitrary points

	void eval_at( const Location& loc, ValueType& res ) const {
		typename ShapeFuncType::Interpolation itp ;
		m_shape.interpolate( loc, itp );

		d6::set_zero( res ) ;
		ValueType node_val ;
		for( Index k = 0 ; k < itp.nodes.rows() ; ++k ) {
			eval_at_node< ValueType >(  itp.nodes[k], Segmenter<D, ValueType>::val2seg( node_val ) ) ;
			res += itp.coeffs[k] * node_val ;
		}
	}

	ValueType eval_at( const Location& loc ) const {
		ValueType seg ; eval_at( loc, seg );
		return seg ;
	}

	ValueType operator() ( const Location& loc ) const { return eval_at(loc) ; }
	template < typename SFINAEShapeFunc = ShapeFuncType >
	typename std::enable_if<SFINAEShapeFunc::is_mesh_based, ValueType>::type
	operator() ( const Vec& x ) const { return eval_at(m_shape.locate(x)) ; }

	template< typename DestShape >
	typename std::enable_if< std::is_same<typename DestShape::DOFDefinition, typename ShapeFuncType::DOFDefinition >::value,
	Interpolation< Derived, D, ShapeFuncT, DestShape > >::type
	interpolate() const
	{ return Interpolation< Derived, D, ShapeFuncT, DestShape >(*this, m_shape.dofDefinition() ) ; }

	template< typename DestShape >
	Interpolation< Derived, D, ShapeFuncT, DestShape > interpolate( const typename DestShape::DOFDefinition& dest ) const
	{ return Interpolation< Derived, D, ShapeFuncT, DestShape >(*this, dest) ; }

	// Abstract

	Derived& derived()
	{ return static_cast< Derived& >( *this ) ; }
	const Derived& derived() const
	{ return static_cast< const Derived& >( *this ) ; }

	Index size( ) const
	{ return derived().size() ; }

	template < typename Agg >
	void eval_at_node( Index i, typename Segmenter<D, Agg>::Seg v ) const
	{ derived().template eval_at_node< Agg >(i, v) ; }

	const ShapeFuncType& shape() const { return m_shape ; }

protected:
	ShapeFuncT m_shape ;
};

template < typename Derived, Index D, typename ShapeFuncT, typename DestShapeFuncT >
struct NonTrivialInterpolation {
	typedef FieldFuncBase< Derived, D, ShapeFuncT > Func ;

	const Func& src ;
	const typename DestShapeFuncT::DOFDefinition& dest ;

	explicit NonTrivialInterpolation( const Func& func_, const typename DestShapeFuncT::DOFDefinition& dest_ )
		: src( func_ ), dest(dest_) {}
};

template < typename Derived, Index D, typename ShapeFuncT, typename DestShapeFuncT >
struct Interpolation : public NonTrivialInterpolation< Derived, D, ShapeFuncT, DestShapeFuncT >  {
	typedef NonTrivialInterpolation< Derived, D, ShapeFuncT, DestShapeFuncT > Base ;
	explicit Interpolation( const typename Base::Func& func, const typename DestShapeFuncT::DOFDefinition& dest )
		:	Base( func, dest )
	{}
};


// Trivial specializationfor non-interpolating interpolation
/* Nope DofDefinition may de different
template < typename Derived, Index D, typename ShapeFuncT >
struct Interpolation< Derived, D, ShapeFuncT, ShapeFuncT >
		: public FieldFuncBase< Interpolation<Derived, D, ShapeFuncT, ShapeFuncT >, D, ShapeFuncT >
{
	typedef FieldFuncBase< Interpolation, D, ShapeFuncT > Base ;
	typedef FieldFuncBase< Derived      , D, ShapeFuncT > Func ;
	const Func& func ;

	explicit Interpolation( const Func& func_, const typename ShapeFuncT::DOFDefinition& )
		:  Base( func_.shape() ), func( func_ )
	{}

	template < typename Agg >
	void eval_at_node( Index i, typename Segmenter<D, Agg>::Seg v ) const
	{
		return func.eval_at_node< Agg >( i, v );
	}
	Index size( ) const
	{
		return func.size() ;
	}
}; */


} //d6


#endif
