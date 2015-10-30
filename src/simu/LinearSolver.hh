#ifndef D6_LINEAR_SOLVER_HH
#define D6_LINEAR_SOLVER_HH

#include "utils/alg.hh"

#include "FormBuilder.hh"
#include <bogus/Core/Block.fwd.hpp>


namespace d6 {

template < typename Derived, typename OtherDerived >
Scalar solveSDP( const bogus::SparseBlockMatrixBase< Derived >& M,
						  const bogus::SparseBlockMatrixBase< OtherDerived >& P,
						  const DynVec &lhs,
						  DynVec &res  ) ;


} //d6


#endif
