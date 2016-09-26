#include "mono/Phase.hh"
#include "diphasic/FluidPhase.hh"

#include "visu/Offline.hh"

#include "utils/Log.hh"
#include "utils/Config.hh"
#include "utils/units.hh"

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

	double convol_width = 3. * offline.box()[1]/(res-1) ;
	int convol_samp = 5 ;

	const Config& c = offline.config() ;

	for( unsigned i = 0 ; i < res ; ++i ) {
		Vec x ;

		Scalar phi = 0, lambda = 0, p = 0 ;
		const Scalar base =  i * ( offline.box()[1]/(res-1)) ;

		for( int k = -convol_samp ; k <= convol_samp ; ++k ) {
			const Scalar convol_x = k * convol_width / convol_samp ;
			const Scalar convol_w = 1./(2*convol_samp + 1 ) ;

			x[1] = base + convol_x ;

			for( unsigned j = 0 ; j < res ; ++j ) {
				x[0] = j * ( offline.box()[0]/(res-1)) ;
				const Vec& xc = offline.meshes().primal().clamp_point( x ) ;
//				const Vec& xc = Vec::Constant(1.e-4).cwiseMax(x).cwiseMin(offline.box() - Vec::Constant(1.e-4) ) ;

				phi    += convol_w * offline.grains().fraction( xc ) ;
				lambda += convol_w * offline.grains().stresses.trace().eval_at( offline.meshes().dual().locate( xc ) ) ;
				p      += convol_w * offline.fluid().pressure( xc ) ;
			}

		}
		phi /= res ;
		lambda /= res ;
		p /= res ;



		const Scalar sL = offline.box()[1] ;
		const Scalar sP = c.units().fromSI(Units::Stress)
		        * ( c.gravity.norm() * c.units().toSI(Units::Acceleration)  )
		        * c.fluidVolMass * c.units().toSI(Units::VolumicMass)
		        * c.box[1] * c.units().toSI(Units::Length);

		out << base / sL  << "\t" << phi << "\t"
		    << lambda / (sP * offline.config().alpha() * offline.config().phiMax ) << "\t"
		    << p / ( sP  )
		    << "\n" ;
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
