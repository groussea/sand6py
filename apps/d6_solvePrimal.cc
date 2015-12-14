
#include "utils/string.hh"
#include "utils/Log.hh"

#include "solve/Primal.hh"
#include "solve/PrimalData.hh"

#include <bogus/Core/Block.impl.hpp>
#include <bogus/Core/BlockSolvers/ProjectedGradient.hpp>

using namespace d6 ;

int main( int argc, char* argv[] ) {

	const char * problem = nullptr ;

	Primal::SolverOptions options ;


	for( int i = 1 ; i < argc ; ++i )
	{
		if( argv[i][0] == '-' ){
			switch(argv[i][1]) {
			case 'n':
				if( ++i == argc ) break ;
				options.maxIterations = to_uint( argv[i] ) ;
				break ;
			case 'f':
				if( ++i == argc ) break ;
				options.maxOuterIterations = to_uint( argv[i] ) ;
				break ;
			case 'a':
				if( ++i == argc ) break ;
				options.algorithm = (Primal::SolverOptions::Algorithm) to_uint( argv[i] ) ;
				break ;
			case 'v':
				if( ++i == argc ) break ;
				options.projectedGradientVariant = to_int( argv[i] ) ;
				break ;
			case 't':
				if( ++i == argc ) break ;
				options.tolerance = to_double( argv[i] ) ;
				break;
			}
		} else {
			problem = argv[i] ;
		}
	}

	Log::Verbose() << "Using algorithm: " ;
	switch( options.algorithm ) {
	case Primal::SolverOptions::GaussSeidel:
		Log::Verbose() << "Gauss-Seidel " ;
		break ;
	case Primal::SolverOptions::Cadoux_GS:
		Log::Verbose() << "Cadoux / Gauss-Seidel " ;
		break ;
	case Primal::SolverOptions::Cadoux_PG_NoAssembly:
		Log::Verbose() << "[No Assembly] " ;
	default:
		Log::Verbose() << "Cadoux / Projected Gradient -- " ;

		switch( (bogus::projected_gradient::Variant) options.projectedGradientVariant ) {
		case bogus::projected_gradient::Descent:
			Log::Verbose() << "Descent" ;
			break;
		case bogus::projected_gradient::Conjugated:
			Log::Verbose() << "Conjugated" ;
			break;
		case bogus::projected_gradient::APGD:
			Log::Verbose() << "APGD" ;
			break;
		case bogus::projected_gradient::SPG:
			Log::Verbose() << "SPG" ;
			break;
		default:
			Log::Verbose() << "Standard" ;
			break;
		}

		break ;

	}
	Log::Verbose() << std::endl ;


	PrimalData data ;
	if( data.load( problem ) ) {

		Primal::SolverStats stats ;

		DynVec r ( data.n() * 6 ) ;
		r.setZero() ;

		Primal solver( data ) ;
		solver.solve( options, r, stats ) ;

		Log::Info() << arg3( "Primal: %1 iterations,\t err= %2,\t time= %3 ",
							 stats.nIterations, stats.residual, stats.time ) << std::endl ;

		return 0 ;
	}

	return 1 ;
}
