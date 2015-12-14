#ifndef D6_PRIMAL_HH
#define D6_PRIMAL_HH

#include "utils/alg.hh"

namespace  d6 {

struct PrimalData ;

class Primal {
public:

	struct SolverStats
	{
		Scalar residual ;
		Scalar time ;
		unsigned nIterations ;

		// For internal use
		void ackResidual( unsigned iter, Scalar err ) ;
	} ;

	struct SolverOptions
	{
		enum Algorithm {
			GaussSeidel,
			Cadoux_GS,
			Cadoux_PG,
			Cadoux_PG_NoAssembly
		};

		Algorithm algorithm ;

		unsigned maxIterations ;	  //Inner
		unsigned maxOuterIterations ; //For Cadoux algorithm

		int projectedGradientVariant ;

		Scalar tolerance ;

		SolverOptions() ;
	};

	Primal( const PrimalData &data ) ;

	Scalar solve( const SolverOptions &options, DynVec& lambda, SolverStats &stats ) const ;

private:

	const PrimalData& m_data ;

};

} //d6

#endif
