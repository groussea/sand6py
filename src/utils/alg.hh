#ifndef D6_ALG_HH
#define D6_ALG_HH

#include "scalar.hh"

#include <Eigen/Core>

namespace d6 {

typedef Eigen::Matrix< Scalar, Eigen::Dynamic, 1 > DynVec;
typedef Eigen::Array < Scalar, Eigen::Dynamic, 1 > DynArr;
typedef Eigen::Matrix< Scalar, Eigen::Dynamic, Eigen::Dynamic > DynMat;

typedef Eigen::Matrix< Scalar, 3, 3 > Mat ;
typedef Eigen::Matrix< Scalar, 3, 1 > Vec ;
typedef Eigen::Array < Scalar, 3, 1 > Arr ;

typedef Eigen::Matrix< Index, 3, 1 > Vec3i ;
typedef Eigen::Array < Index, 3, 1 > Arr3i ;

typedef Eigen::Matrix< Scalar, 6, 6 > Mat66 ;
typedef Eigen::Matrix< Scalar, 6, 1 > Vec6 ;

typedef Eigen::Matrix< Scalar, 3, Eigen::Dynamic > DynMat3 ;
typedef Eigen::Matrix< Scalar, 6, Eigen::Dynamic > DynMat6 ;

template < Index Dimension >
struct Segmenter {
	typedef Eigen::Matrix< Scalar, Dimension, 1 > ValueType ;

	typedef typename DynVec::template FixedSegmentReturnType< Dimension >::Type Seg ;
	typedef typename DynVec::template ConstFixedSegmentReturnType< Dimension >::Type ConstSeg ;

	static inline Seg segment( DynVec& vec, const Index i ) {
		return vec.template segment< Dimension >( i*Dimension ) ;
	}
	static inline ConstSeg segment( const DynVec& vec, const Index i ) {
		return vec.template segment< Dimension >( i*Dimension ) ;
	}

	Seg val2seg( ValueType & v) { return v.segment< Dimension >(0) ; }

};
template < >
struct Segmenter< 1 > {
	typedef Scalar ValueType ;

	typedef Scalar& Seg ;
	typedef const Scalar& ConstSeg ;

	static inline Seg segment( DynVec& vec, const Index i ) {
		return vec[ i ] ;
	}
	static inline ConstSeg segment( const DynVec& vec, const Index i ) {
		return vec[ i ] ;
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
