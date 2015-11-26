#ifndef D6_FIELD_BASE_HH
#define D6_FIELD_BASE_HH

#include "Expr.hh"
#include "FieldFuncBase.hh"

namespace d6 {

template< typename Derived >
struct FieldTraits
{

};


template< typename Derived >
class FieldBase : public Expr< typename FieldTraits< Derived >::ValueType >
{

public:

	typedef FieldTraits< Derived > Traits ;
	typedef MeshBase< typename Traits::MeshType > MeshType ;

	static constexpr Index D = Traits::Dimension ;
	typedef typename Traits::ValueType ValueType ;

	typedef typename Segmenter< D >::Seg Seg  ;
	typedef typename Segmenter< D >::ConstSeg ConstSeg  ;

	explicit FieldBase( const MeshType& mesh )
		: m_mesh( mesh ), m_size( mesh.nNodes() )
	{
		m_data.resize( D * m_size );
	}

	const DynVec& flatten() const { return m_data ; }
	DynVec& flatten() { return m_data ; }

	const MeshType& mesh() const { return m_mesh ; }

	//Global setters

	Derived& set_zero() ;
	Derived& set_constant( const ValueType& val ) ;

	template < typename Func >
	Derived& operator= ( const FieldFuncBase< Func, D, typename Traits::MeshType > & f )
	{
		assert( f.size() == size() ) ;
		#pragma omp parallel for
		for( Index  i = 0 ; i < size() ; ++i ) {
			f.eval_at_node( i, segment(i) );
		}
		return derived();
	}

	// Interpolation

	void  add_at( const Vec& x, const ValueType& val ) ;
	void  add_at( const typename MeshType::Interpolation &itp, const ValueType& val ) ;
	void eval_at( const Vec& x, ValueType& res ) const ;

	ValueType eval_at( const Vec& x ) const {
		ValueType seg ; eval_at( x, seg );
		return seg ;
	}

	// Value at node
	Seg      segment( const Index i ) { return Segmenter< D >::segment( m_data, i) ; }
	ConstSeg segment( const Index i ) const { return Segmenter< D >::segment(m_data, i) ; }

	Index size() const { return m_size ; }

	// Operators
	ValueType operator() ( const Vec&  x ) const { return eval_at(x) ; }
	ConstSeg  operator[] ( const Index i ) const { return segment(i) ; }
	Seg       operator[] ( const Index i )       { return segment(i) ; }

	// Info
	Scalar max_abs() const { return m_data.lpNorm< Eigen::Infinity >() ; }

	//Operations

	Derived& multiply_by( const ScalarField& field ) ;
	Derived&   divide_by( const ScalarField& field ) ;
	Derived&   divide_by_positive( const ScalarField& field, Scalar min = 1.e-16 ) ;

	// Serialization

	template < typename Archive >
	void serialize( Archive &ar, unsigned int ) {
		ar & m_size ;
		ar & m_data ;
	}

	Derived& derived()
	{ return static_cast< Derived& >( *this ) ; }
	const Derived& derived() const
	{ return static_cast< const Derived& >( *this ) ; }

protected:
	const MeshType & m_mesh ;
	Index m_size ;
	DynVec   m_data ;

} ;

} //d6

#endif

