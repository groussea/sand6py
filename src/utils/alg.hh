#ifndef D6_ALG_HH
#define D6_ALG_HH

#include "scalar.hh"

#include <Eigen/Core>

namespace d6 {

typedef Eigen::Matrix< Scalar, Eigen::Dynamic, 1 > DynVec;
typedef Eigen::Matrix< Scalar, Eigen::Dynamic, Eigen::Dynamic > DynMat;

typedef Eigen::Matrix< Scalar, 3, 3 > Mat ;
typedef Eigen::Matrix< Scalar, 3, 1 > Vec ;

typedef Eigen::Matrix< Index, 3, 1 > Vec3i ;

typedef Eigen::Matrix< Scalar, 6, 6 > Mat66 ;
typedef Eigen::Matrix< Scalar, 6, 1 > Vec6 ;

typedef Eigen::Matrix< Scalar, 3, Eigen::Dynamic > DynMat3 ;
typedef Eigen::Matrix< Scalar, 6, Eigen::Dynamic > DynMat6 ;

template < Index Dimension >
struct Segmenter {

	typedef typename DynVec::template FixedSegmentReturnType< Dimension >::Type Seg ;
	typedef typename DynVec::template ConstFixedSegmentReturnType< Dimension >::Type ConstSeg ;

	static inline Seg segment( DynVec& vec, const Index i ) {
		return vec.template segment< Dimension >( i*Dimension ) ;
	}
	static inline ConstSeg segment( const DynVec& vec, const Index i ) {
		return vec.template segment< Dimension >( i*Dimension ) ;
	}
};
template < >
struct Segmenter< 1 > {
	typedef Scalar& Seg ;
	typedef const Scalar& ConstSeg ;

	static inline Seg segment( DynVec& vec, const Index i ) {
		return vec[ i ] ;
	}
	static inline ConstSeg segment( const DynVec& vec, const Index i ) {
		return vec[ i ] ;
	}

} ;

template < typename Derived >
inline void set_zero( Eigen::MatrixBase< Derived > &mat ) {
	mat.setZero() ;
}
inline void set_zero( Scalar &s ) {
	s = 0 ;
}

} //d6

#endif
