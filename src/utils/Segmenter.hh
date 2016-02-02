#ifndef D6_SEGMENTER_HH
#define D6_SEGMENTER_HH

#include "scalar.hh"

#include <Eigen/Core>

namespace d6 {

typedef Eigen::Matrix< Scalar, Eigen::Dynamic, 1 > DynVec;
typedef Eigen::Array < Scalar, Eigen::Dynamic, 1 > DynArr;
typedef Eigen::Matrix< Scalar, Eigen::Dynamic, Eigen::Dynamic > DynMat;

template < Index Dimension, typename Aggregate = DynVec >
struct Segmenter {
	typedef Eigen::Matrix< Scalar, Dimension, 1 > ValueType ;

	typedef typename Aggregate::template FixedSegmentReturnType< Dimension >::Type Seg ;
	typedef typename Aggregate::template ConstFixedSegmentReturnType< Dimension >::Type ConstSeg ;

	static inline Seg segment( Aggregate& vec, const Index i ) {
		return vec.template segment< Dimension >( i*Dimension ) ;
	}
	static inline ConstSeg segment( const Aggregate& vec, const Index i ) {
		return vec.template segment< Dimension >( i*Dimension ) ;
	}
	static inline Seg head( Aggregate& vec ) {
		return vec.template head< Dimension >( ) ;
	}
	static inline ConstSeg head( const Aggregate& vec ) {
		return vec.template head< Dimension >( ) ;
	}
	static inline Seg tail( Aggregate& vec ) {
		return vec.template tail< Dimension >( ) ;
	}
	static inline ConstSeg tail( const Aggregate& vec ) {
		return vec.template tail< Dimension >( ) ;
	}

	Seg val2seg( ValueType & v) { return v.template head< Dimension >() ; }

};
template < typename Aggregate >
struct Segmenter< 1, Aggregate > {
	typedef Scalar ValueType ;

	typedef Scalar& Seg ;
	typedef const Scalar& ConstSeg ;

	static inline Seg segment( Aggregate& vec, const Index i ) {
		return vec[ i ] ;
	}
	static inline ConstSeg segment( const Aggregate& vec, const Index i ) {
		return vec[ i ] ;
	}
	static inline Seg head( Aggregate& vec ) {
		return vec[0] ;
	}
	static inline ConstSeg head( const Aggregate& vec ) {
		return vec[0] ;
	}
	static inline Seg tail( Aggregate& vec ) {
		return vec[vec.rows()-1] ;
	}
	static inline ConstSeg tail( const Aggregate& vec ) {
		return vec[vec.rows()-1] ;
	}

	Seg val2seg( ValueType& v) { return v ; }
} ;


template < typename Derived >
inline void set_zero( Eigen::MatrixBase< Derived > &mat ) {
	mat.setZero() ;
}
inline void set_zero( Scalar &s ) {
	s = 0 ;
}


template < Index Dimension, typename Derived >
void mul_compwise( DynVec& vec, const Eigen::MatrixBase<Derived> & scalar ) {
	assert( vec.rows() == scalar.rows() * Dimension ) ;
	Eigen::Matrix< Scalar, Dimension, Eigen::Dynamic >::Map( vec.data(), Dimension, scalar.rows() )
			*= scalar.asDiagonal() ;
}
template < Index Dimension, typename Derived >
void div_compwise( DynVec& vec,  const Eigen::MatrixBase<Derived> & scalar ) {
	assert( vec.rows() == scalar.rows() * Dimension ) ;
	Eigen::Matrix< Scalar, Dimension, Eigen::Dynamic >::Map( vec.data(), Dimension, scalar.rows() )
			*= (1./scalar.array()).matrix().asDiagonal() ;
}
template < Index Dimension, typename Derived >
void set_compwise( DynVec& vec, const Eigen::MatrixBase< Derived > & scalar ) {
	vec.setOnes( Dimension * scalar.rows() ) ;
	mul_compwise< Dimension >( vec, scalar ) ;
}

template < Index Dimension >
typename Eigen::Matrix< Scalar, Dimension, Eigen::Dynamic >::MapType::RowXpr component( DynVec& vec, Index row ) {
	 return Eigen::Matrix< Scalar, Dimension, Eigen::Dynamic >::Map( vec.data(), Dimension, vec.rows()/Dimension ).row( row ) ;
 }
template < Index Dimension >
typename Eigen::Matrix< Scalar, Dimension, Eigen::Dynamic >::ConstMapType::ConstRowXpr component( const DynVec& vec, Index row ) {
	 return Eigen::Matrix< Scalar, Dimension, Eigen::Dynamic >::Map( vec.data(), Dimension, vec.rows()/Dimension ).row( row ) ;
 }

} //d6

#endif
