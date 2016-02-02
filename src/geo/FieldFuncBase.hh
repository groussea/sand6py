#ifndef D6_FIELD_FUNC_HH
#define D6_FIELD_FUNC_HH

#include "utils/alg.hh"
#include "geo.fwd.hh"

namespace d6 {

template < typename Derived, Index D, typename ShapeFuncT >
struct FieldFuncBase
{
	typedef ShapeFuncBase< ShapeFuncT > ShapeFuncType ;

	typedef typename Segmenter< D >::ValueType ValueType ;
	typedef typename Segmenter< D >::Seg Seg  ;

	explicit FieldFuncBase( const ShapeFuncType& shape ) : m_shape(shape.derived()) {}


	ValueType operator[]( Index i ) const
	{
		ValueType v ;
		return eval_at_node( i, Segmenter< D >::val2seg(v) ) ;
	}

	void eval_at_node( Index i, Seg v ) const
	{
		return static_cast< const Derived& >(*this).eval_at_node(i, v) ;
	}
	Index size( ) const
	{
		return static_cast< const Derived& >(*this).size() ;
	}

	const ShapeFuncType& shape() const { return m_shape ; }

protected:
	ShapeFuncT m_shape ;
};


} //d6


#endif
