#ifndef D6_BLOCK_MAT_HH
#define D6_BLOCK_MAT_HH

#include "utils/alg.hh"
#include <bogus/Core/Block.hpp>

namespace d6 {

template < Index Rows, Index Cols >
struct FormMat {
	typedef Eigen::Matrix< Scalar, Rows, Cols > BlockT ;
	typedef bogus::SparseBlockMatrix< BlockT > Type ;
	typedef bogus::SparseBlockMatrix< BlockT, bogus::SYMMETRIC > SymType ;
};

} //d6

#endif
