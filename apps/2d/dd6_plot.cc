#include "mono/Phase.hh"
#include "diphasic/FluidPhase.hh"

#include "visu/Offline.hh"

#include "utils/Log.hh"

#include <iostream>

using namespace d6 ;

static void usage( const char *name )
{
	std::cout << "Usage: " << name
	          << " [sim_dir=out] [options] "
	          << "\nGenerates average plots."
	          << "\n\n" ;

	std::cout << "Options:\n"
	          << "-? \t Display this help message and exit\n"
	          << "-n frame_id \t Stop at frame frame_id\n"
	          << std::endl ;
}

bool plot( unsigned frame, unsigned res, Offline& offline, std::ostream& out ) {

	if(!offline.load_frame(frame))
		return false ;

	for( unsigned i = 0 ; i < res ; ++i ) {
		Vec x ;
		x[1] = i * ( offline.box()[1]/(res-1)) ;

		Scalar phi = 0, lambda = 0, p = 0 ;

		for( unsigned j = 0 ; j < res ; ++j ) {
			x[0] = j * ( offline.box()[0]/(res-1)) ;

			phi    += offline.grains().fraction( x ) ;
			lambda += offline.grains().stresses.trace().eval_at( offline.meshes().dual().locate( x ) ) ;
			p      += offline.fluid().pressure( x ) ;
		}
		phi /= res ;
		lambda /= res ;
		p /= res ;

		out << x[1] << "\t" << phi << "\t"
		                  << lambda<< "\t" << p << "\n" ;
	}

	return true ;
}

int main( int argc, const char* argv[] ) {


	const char * base_dir = "out" ;
	unsigned last_frame = -1 ;
	unsigned res = 100 ;

	for( int i = 1 ; i < argc ; ++i )
	{
		if( argv[i][0] == '-' ){
			switch(argv[i][1]) {
			case '?':
				usage(argv[0]) ;
				return 0 ;
			case 'n':
				if( ++i == argc ) break ;
				last_frame = d6::to_uint( argv[i] ) ;
				break ;
			case 'r':
				if( ++i == argc ) break ;
				res = d6::to_uint( argv[i] ) ;
				break ;
			}
		} else {
			base_dir = argv[i] ;
		}
	}

	Log::Config::get().level = Log::L_Error ; ;

	Offline offline( base_dir ) ;
//	const Config &config = offline.config() ;

//	for( unsigned cur_frame = 0 ;
//	     ( last_frame > cur_frame ) && offline.load_frame( cur_frame ) ;
//	     ++cur_frame )
//	{
//		const Scalar t =config.time( cur_frame ) ;
//		const Scalar vol = remaining_volume( offline.particles(), config.box )  ;
//		std::cout << cur_frame << "\t" << t << "\t" << (t * config.units().toSI( Units::Time ))
//		          << "\t" << vol << "\t" << (vol * config.units().toSI( Units::Volume ) )
//		          << "\n" ;
//	}
	if( plot( last_frame, res, offline, std::cout ) )
		return 0 ;

	return 1 ;

}
