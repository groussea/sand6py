#ifndef D6_FIELD_FUNC_HH
#define D6_FIELD_FUNC_HH

#include "utils/alg.hh"
#include "geo.fwd.hh"

namespace d6 {

template < typename Derived, Index D, typename MeshT >
struct FieldFuncBase
{
	typedef MeshBase< MeshT > MeshType ;

	typedef typename Segmenter< D >::ValueType ValueType ;
	typedef typename Segmenter< D >::Seg Seg  ;

	explicit FieldFuncBase( const MeshType& mesh ) : m_mesh(mesh) {}


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

	const MeshType& mesh() const { return m_mesh ; }

protected:
	const MeshType& m_mesh ;
};


} //d6


#endif
