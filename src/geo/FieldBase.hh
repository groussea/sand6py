#ifndef D6_FIELD_BASE_HH
#define D6_FIELD_BASE_HH

#include "FieldFuncBase.hh"
#include "ShapeFunctionBase.hh"

namespace d6 {

template< typename Derived >
struct FieldTraits
{

};



template< typename Derived >
class FieldBase : public FieldFuncBase< Derived,
		FieldTraits< Derived >::Dimension,
		typename FieldTraits< Derived >::ShapeFuncImpl >
{

public:

	typedef FieldTraits< Derived > Traits ;

	typedef typename Traits::ShapeFuncImpl ShapeFuncImpl ;
	typedef ShapeFuncBase< ShapeFuncImpl > ShapeFuncType ;
	typedef typename ShapeFuncType::Location Location ;

	static constexpr Index D = Traits::Dimension ;

	typedef FieldFuncBase< Derived, D, ShapeFuncImpl > Base ;
	using Base::m_shape ;
	using Base::derived ;

	typedef typename Segmenter< D >::ValueType ValueType ;
	typedef typename Segmenter< D >::Seg Seg  ;
	typedef typename Segmenter< D >::ConstSeg ConstSeg  ;

	typedef AbstractScalarField< ShapeFuncImpl > ScalarField ;

	explicit FieldBase( const ShapeFuncType& shape )
		: Base( shape )
	{
		fit_shape() ;
	}

	void fit_shape() {
		m_size = m_shape.nDOF() ;
		m_data.resize( D * m_size );
	}

	const DynVec& flatten() const { return m_data ; }
	DynVec& flatten() { return m_data ; }

	const ShapeFuncImpl& shape() const { return m_shape ; }

	//Global setters

	Derived& set_zero() ;
	Derived& set_constant( const ValueType& val ) ;

	template < typename Func >
	Derived& operator= ( const FieldFuncBase< Func, D, ShapeFuncImpl > & f )
	{
		assert( f.size() == size() ) ;
		#pragma omp parallel for
		for( Index  i = 0 ; i < size() ; ++i ) {
			f.template eval_at_node< DynVec >( i, segment(i) );
		}
		return derived();
	}

	// Accumulation, interpolation

	void  add_at( const Location& loc, const ValueType& val ) ;
	void  add_at( const typename ShapeFuncType::Interpolation &itp, const ValueType& val ) ;

	template < typename Func, typename OtherShape >
	Derived& from_interpolation( const FieldFuncBase< Func, D, OtherShape > &f ) ;

	// Value at node
	Seg      segment( const Index i ) { return Segmenter< D >::segment( m_data, i) ; }
	ConstSeg segment( const Index i ) const { return Segmenter< D >::segment(m_data, i) ; }

	template < typename Agg >
	void eval_at_node( Index i, typename Segmenter<D, Agg>::Seg v ) const {	v = segment(i) ; }

	Index size() const { return m_size ; }

	// Operators

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


protected:

	template < typename Func, typename OtherShape >
	typename std::enable_if< OtherShape::is_mesh_based>::type
	integrate( const FieldFuncBase< Func, D, OtherShape > &f ) ;

	template < typename Func, typename OtherShape >
	typename std::enable_if<!OtherShape::is_mesh_based >::type
	integrate( const FieldFuncBase< Func, D, OtherShape > &f ) ;

	Index m_size ;
	DynVec   m_data ;

} ;

#define D6_MAKE_FIELD_CONSTRUCTORS_AND_ASSIGNMENT_OPERATORS( FieldName ) \
	explicit FieldName( const typename Base::ShapeFuncType& shape ) \
	: Base( shape ) \
	{} \
	explicit FieldName( const typename Base::ShapeFuncType::DOFDefinition& mesh ) \
	: Base( typename Base::ShapeFuncImpl( mesh ) ) \
	{} \
	/*! Assignment from field func */\
	template <typename Func> \
	FieldName( const FieldFuncBase< Func, Base::D, typename Base::ShapeFuncImpl > & func ) \
		: Base( func.shape() ) \
	{ Base::operator=( func ); } \
	/*! Constructor from field func */\
	template <typename Func> \
	FieldName& operator= ( const FieldFuncBase< Func, Base::D, typename Base::ShapeFuncImpl> & func ) \
	{ return Base::operator=( func ); } \
	/*! Assignment from non-trivial interpolation */\
	template < typename Func, typename OtherShape > \
	FieldName& operator= ( const NonTrivialInterpolation< Func, Base::D, OtherShape, typename Base::ShapeFuncImpl > & f ) \
	{ return Base::from_interpolation ( f.func ) ; } \
	/*! Constructors from non-trivial interpolation */\
	template < typename Func, typename OtherShape > \
	FieldName ( const NonTrivialInterpolation< Func, Base::D, OtherShape, typename Base::ShapeFuncImpl > & f ) \
	: Base(  typename Base::ShapeFuncImpl( f.func.shape().dofDefinition() ) ) \
	{ Base::from_interpolation ( f.func ) ; } \
	template < typename Func, typename OtherShape > \
	FieldName ( const typename Base::ShapeFuncType::DOFDefinition& mesh, \
				const NonTrivialInterpolation< Func, Base::D, OtherShape, typename Base::ShapeFuncImpl > & f ) \
	: Base( typename Base::ShapeFuncImpl( mesh ) ) \
	{ Base::from_interpolation ( f.func ) ; } \


} //d6

#endif

