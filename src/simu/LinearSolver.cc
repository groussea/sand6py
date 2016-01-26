#include "LinearSolver.hh"

#include "FormBuilder.hh"

#include "utils/Log.hh"

#include <bogus/Core/Block.impl.hpp>
#include <bogus/Core/BlockSolvers.impl.hpp>

//#define USE_EIGEN
#ifdef USE_EIGEN
#include <bogus/Core/Block.io.hpp>
#include <Eigen/SparseCholesky>
#include <Eigen/Cholesky>
#endif

namespace d6 {

template < typename Derived, typename OtherDerived >
Scalar solveSDP( const bogus::SparseBlockMatrixBase< Derived >& M,
						  const bogus::SparseBlockMatrixBase< OtherDerived >& P,
						  const DynVec &rhs,
						  DynVec &res  )
{

#ifdef USE_EIGEN
	(void) P ;

	typedef Eigen::SparseMatrix< Scalar > SM ;
	SM csr ;
	bogus::convert( M, csr ) ;

	csr.prune(1.) ;

	Eigen::SimplicialLDLT< SM > ldlt( csr ) ;
	std::cerr << "DLLD " << ldlt.info() << std::endl ;
	res = ldlt.solve( rhs ) ;

	return 0 ;
#else
	typedef bogus::MatrixPreconditioner<  OtherDerived                   > Precond ;

	Scalar cgres = 0. ;

	bogus::Krylov< Derived,	Precond::template Type > cg ( M.derived() ) ;
	cg.preconditioner().setPreconditionerMatrix( P.derived() ) ;

	cgres = cg.asCG().solve( rhs, res ) ;
	Log::Debug() << "CG res: " << cgres << std::endl ;


	return cgres ;
#endif
}

template Scalar solveSDP( const bogus::SparseBlockMatrixBase< typename FormMat<WD,WD>::Type >& ,
const bogus::SparseBlockMatrixBase< typename FormMat<WD,WD>::SymType >&,
const DynVec &,
DynVec &) ;

} //d6
