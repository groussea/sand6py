#ifndef D6_FIELD_BASE_HH
#define D6_FIELD_BASE_HH

#include "utils/alg.hh"
#include "geo.fwd.hh"

namespace d6 {

template< typename Derived >
struct FieldTraits
{

};

template< typename Derived >
struct FieldBase
{

public:

	typedef FieldTraits< Derived > Traits ;
	typedef MeshBase< typename Traits::MeshType > MeshType ;

	static constexpr Index D = Traits::Dimension ;
	typedef typename Traits::ValueType ValueType ;

	typedef typename Segmenter< D >::Seg Seg  ;
	typedef typename Segmenter< D >::ConstSeg ConstSeg  ;


	FieldBase( const MeshType& mesh )
		: m_mesh( mesh ), m_size( mesh.nNodes() )
	{
		m_data.resize( D * m_size );
	}

	const DynVec& flatten() const { return m_data ; }
	DynVec& flatten() { return m_data ; }

	//
	void set_zero() ;
	void set_constant( const ValueType& val ) ;

	// Interpolation

	void  add_at( const Vec& x, const ValueType& val ) ;
	void eval_at( const Vec& x, ValueType& res ) const ;

	ValueType eval_at( const Vec& x ) const {
		ValueType seg ; eval_at( x, seg );
		return seg ;
	}

	// Value at node
	Seg      segment( const Index i ) { return Segmenter< D >::segment( m_data, i) ; }
	ConstSeg segment( const Index i ) const { return Segmenter< D >::segment(m_data, i) ; }

	Index size() { return m_size ; }

	// Operators
	ValueType operator() ( const Vec&  x ) const { return eval_at(x) ; }
	ConstSeg  operator[] ( const Index i ) const { return segment(i) ; }
	Seg       operator[] ( const Index i )       { return segment(i) ; }

protected:
	const MeshType & m_mesh ;
	const Index m_size ;
	DynVec   m_data ;

} ;

} //d6

#endif

