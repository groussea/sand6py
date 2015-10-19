#ifndef D6_ALG_HH
#define D6_ALG_HH

#include "scalar.hh"

#include <Eigen/Core>

namespace d6 {

typedef Eigen::Matrix< Scalar, Eigen::Dynamic, 1 > DynVec;
typedef Eigen::Matrix< Scalar, Eigen::Dynamic, Eigen::Dynamic > DynMat;

typedef Eigen::Matrix< Scalar, 3, 3 > Mat ;
typedef Eigen::Matrix< Scalar, 3, 1 > Vec ;

typedef Eigen::Matrix< int, 3, 1 > Vec3i ;

typedef Eigen::Matrix< Scalar, 6, 6 > Mat66 ;
typedef Eigen::Matrix< Scalar, 6, 1 > Vec6 ;

}

#endif
