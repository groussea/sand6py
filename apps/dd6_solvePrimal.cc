
#include "utils/string.hh"
#include "utils/Log.hh"

#include "diphasic/DiphasicFriction.hh"

#include <bogus/Core/Block.impl.hpp>
#include <bogus/Core/BlockSolvers/ProjectedGradient.hpp>

using namespace d6 ;

static void usage( const char *name )
{
	std::cout << "Usage: " << name
			  << " problem_file [options] "
			  << "\nOffline DCFP solver. "
			  << "\n\n" ;

	std::cout << "Options:\n"
			  << "-? \t Display this help message and exit\n"
			  << "-a algo \t Algorithm id in [0,3]  \n"
			  << "-g variant \t PG variant id in [0,4]  \n"
			  << "-c useCadoux \t Whether to use Cadoux fp alg \n"
			  << "-t tol \t Tolerance \n"
			  << "-n nIters \t Max number of inner solver iterations  \n"
			  << "-f nIters \t Max number of outer fixed-point iterations  \n"
			  << "-o  \t Use unbiased residual evaluation function (Alart-Curnier)\n"
			  << "-e time \t Abort after given computation time budget\n"
			  << "-v level \t Verbosity level\n"
			  << std::endl ;
}


int main( int argc, char* argv[] ) {

	const char * problem = nullptr ;

	FrictionSolver::Stats stats ;
	DiphasicFrictionSolver::Options options ;

	for( int i = 1 ; i < argc ; ++i )
	{
		if( argv[i][0] == '-' ){
			switch(argv[i][1]) {
			case '?':
				usage(argv[0]) ;
				return 0 ;
			case 'c':
				if( ++i == argc ) break ;
				options.useCadoux = to_bool( argv[i] ) ;
				break ;
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
				options.algorithm = (DiphasicFrictionSolver::Options::Algorithm) to_uint( argv[i] ) ;
				break ;
			case 'g':
				if( ++i == argc ) break ;
				options.projectedGradientVariant = to_int( argv[i] ) ;
				break ;
			case 't':
				if( ++i == argc ) break ;
				options.tolerance = to_double( argv[i] ) ;
				break;
			case 'o':
				stats.shouldLogAC = true ;
				break ;
			case 'e':
				if( ++i == argc ) break ;
				stats.timeOut = to_double( argv[i] ) ;
			case 'v':
				if( ++i == argc ) break ;
				d6::Log::Config::get().setLevel( argv[i] ) ;
			}
		} else {
			problem = argv[i] ;
		}
	}
	if( !problem ) {
		usage(argv[0]) ;
		return 1 ;
	}

	Log::Verbose() << "Using algorithm: " ;
	if( options.useCadoux)
		Log::Verbose() << "Cadoux / " ;

	switch( options.algorithm ) {
	case DiphasicFrictionSolver::Options::ADMM:
		Log::Verbose() << "ADMM" ;
	case DiphasicFrictionSolver::Options::PG_CG_Red:
		Log::Verbose() << "[CG/Red]" ;
		break ;
	case DiphasicFrictionSolver::Options::PG_CG_Stokes:
		Log::Verbose() << "[CG/Stk]" ;
		break ;
	case DiphasicFrictionSolver::Options::PG_Fac_Red:
		Log::Verbose() << "[Fac/Red]" ;
		break ;
	case DiphasicFrictionSolver::Options::PG_Fac_Stokes:
		Log::Verbose() << "[Fac/Stk]" ;
		break ;
	}

	if(options.algorithm != DiphasicFrictionSolver::Options::ADMM)
	{
		Log::Verbose() << " Projected Gradient -- " ;
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
	}

	Log::Verbose() << std::endl ;


	DiphasicPrimalData data ;
	if( data.load( problem ) ) {

		stats.timeOut *= data.n() ;

		DynVec x ( data.s() ) ;
		DynVec r ( data.n() * SD ) ;
		r.setZero() ; x.setZero() ;

		DiphasicFrictionSolver solver( data ) ;
		try {
			solver.solve( options, x, r, stats ) ;
		} catch ( FrictionSolver::Stats::TimeOutException& ){
			Log::Warning() << "Solver exceeded allowed time" << std::endl ;
		}

		Log::Info() << arg3( "Friction: %1 iterations,\t err= %2,\t time= %3 ",
							 stats.nIterations(), stats.residual(), stats.time() ) << std::endl ;


		for( const FrictionSolver::Stats::Entry& e : stats.log() ) {
			std::cout << e.nIterations << "\t" << e.time << "\t" << e.residual << "\n" ;
		}

//		std::cout << x.head( data.m() ).minCoeff()  << "\n"
//				  << x.head( data.m() ).maxCoeff()  << "\n"
//				  << x.head( data.m() ).norm()  << "\n"
//				  << " --- \n "
//				  << x.segment( data.m(), data.r() ).minCoeff()  << "\n"
//				  << x.segment( data.m(), data.r() ).maxCoeff()  << "\n"
//				  << x.segment( data.m(), data.r() ).norm()  << "\n"
//				  << " --- \n "
//				  << x.segment( data.m()+data.r(), data.p() ).minCoeff()  << "\n"
//				  << x.segment( data.m()+data.r(), data.p() ).maxCoeff()  << "\n"
//				  << x.segment( data.m()+data.r(), data.p() ).norm()  << "\n"
//				  << " --- \n "
//				  << r.minCoeff()  << "\n"
//				  << r.maxCoeff()  << "\n"
//				  << r.norm()  << "\n"
//				  << std::endl ;


		return 0 ;
	}

	return 1 ;
}
