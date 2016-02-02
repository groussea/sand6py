#ifndef D6_ALG_HH
#define D6_ALG_HH

#include "Segmenter.hh"

namespace d6 {

#if D6_DIM==2
static constexpr Index WD = 2 ;
static constexpr Index SD = 3 ;
static constexpr Index RD = 1 ;
#else
static constexpr Index WD = 3 ;
static constexpr Index SD = 6 ;
static constexpr Index RD = 3 ;
#endif

typedef Eigen::Matrix< Scalar, WD, WD > Mat ;
typedef Eigen::Matrix< Scalar, WD,  1 > Vec ;
typedef Eigen::Array < Scalar, WD,  1 > Arr ;

typedef Eigen::Matrix< Scalar, SD, SD > MatS ;
typedef Eigen::Matrix< Scalar, SD,  1 > VecS ;
typedef Eigen::Array < Scalar, SD,  1 > ArrS ;

typedef Eigen::Matrix< Scalar, RD, RD > MatR ;
typedef typename Segmenter<RD>::ValueType VecR ;

typedef Eigen::Matrix< Index, WD, 1 > VecWi ;
typedef Eigen::Array < Index, WD, 1 > ArrWi ;

typedef Eigen::Matrix< Scalar, 3, 3 > Mat33 ;
typedef Eigen::Matrix< Scalar, 3, 1 > Vec3 ;
typedef Eigen::Array < Scalar, 3, 1 > Arr3 ;

typedef Eigen::Matrix< Index, 3, 1 > Vec3i ;
typedef Eigen::Array < Index, 3, 1 > Arr3i ;

typedef Eigen::Matrix< Scalar, 6, 6 > Mat66 ;
typedef Eigen::Matrix< Scalar, 6, 1 > Vec6 ;

typedef Eigen::Matrix< Scalar, WD, Eigen::Dynamic > DynMatW ;
typedef Eigen::Matrix< Scalar, SD, Eigen::Dynamic > DynMatS ;
typedef Eigen::Matrix< Scalar, 3, Eigen::Dynamic > DynMat3 ;
typedef Eigen::Matrix< Scalar, 6, Eigen::Dynamic > DynMat6 ;


} //d6

#endif
