#ifndef D6_PRIMAL_HH
#define D6_PRIMAL_HH

#include "utils/alg.hh"
#include <vector>

namespace  d6 {

struct PrimalData ;

class Primal {
public:

	struct SolverStats
	{
		struct Entry {
			Scalar   residual ;
			Scalar   time ;
			unsigned nIterations ;
		};
		struct TimeOutException {} ;
		typedef std::vector< Entry > Log ;

		SolverStats() : shouldLogAC( false ), timeOut( std::numeric_limits<Scalar>::infinity() )
		{}

		Scalar   residual()    const { return m_log.back().residual ; }
		Scalar   time()        const { return m_log.back().time ; }
		unsigned nIterations() const { return m_log.back().nIterations ; }

		void log( unsigned iter, Scalar err, Scalar time ) ;
		const Log& log( )  { return m_log ; }

		//! If true, should store residual computed using Alart-Curnier function
		bool shouldLogAC ;

		//! Throws time-out exception when solver time exceeds timeOut (seconds)
		Scalar timeOut ;
	private:
		Log m_log ;
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
