#include "visu/Offline.hh"

#include "utils/Log.hh"

#include <iostream>

using namespace d6 ;

static Scalar remaining_volume( const Particles& particles, const Vec& box )
{
	const Scalar zmin = box[2] / 2 ;

	Scalar volume = 0 ;
#pragma omp parallel for reduction(+:volume)
	for( unsigned i = 0 ; i < particles.count() ; ++i ) {
		if( particles.centers()(2,i) > zmin ) {
			volume += particles.volumes()[i] ;
		}
	}

	return volume ;
}

int main( int argc, const char* argv[] ) {


	const char * base_dir = "out" ;

	for( int i = 1 ; i < argc ; ++i )
	{
		if( argv[i][0] == '-' ){
		} else {
			base_dir = argv[i] ;
		}
	}

	Offline offline( base_dir ) ;
	const Config &config = offline.config() ;

	Log::Config::get().level = Log::L_Error ; ;

	for( unsigned cur_frame = 0 ;
		 offline.load_frame( cur_frame ) ;
		 ++cur_frame )
	{
		const Scalar t =config.time( cur_frame, 0 ) ;
		const Scalar vol = remaining_volume( offline.particles(), config.box )  ;
		std::cout << cur_frame << "\t" << t << "\t" << (t * config.units().toSI( Units::Time ))
				  << "\t" << vol << "\t" << (vol * config.units().toSI( Units::Volume ) )
				  << "\n" ;
	}

	return 1 ;

}
