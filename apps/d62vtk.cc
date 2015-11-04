
#include "utils/string.hh"

#include "visu/Offline.hh"
#include "visu/VTKParticlesWriter.hh"
#include "visu/VTKFieldWriter.hh"

#include "simu/Phase.hh"

#include "geo/ScalarField.hh"
#include "geo/TensorField.hh"

void dump_frame( const d6::Offline& offline, bool particles,
				 const char* base_dir, unsigned frame )
{

	if(particles) {
		d6::VTKParticlesWriter particlesWriter( base_dir, offline.particles() ) ;
		particlesWriter.startFile( "particles", frame ) ;
		particlesWriter.dump_all() ;
	}

	d6::VTKFieldWriter fieldWriter( base_dir, offline.mesh() ) ;
//	fieldWriter.setMode( d6::VTKWriter::Ascii );
	fieldWriter.startFile( "fields", frame ) ;
	fieldWriter.dump(    "phi", offline.grains().fraction ) ;
	fieldWriter.dump(      "u", offline.grains().velocity ) ;

	d6::ScalarField p   = offline.grains().stresses.trace() ;
	d6::ScalarField dh  = offline.grains().sym_grad.trace() ;
	d6::ScalarField taun= d6::TensorField( offline.grains().stresses.deviatoricPart() ).norm() ;

	fieldWriter.dump(    "p", p ) ;
	fieldWriter.dump(   "dh", dh ) ;
	fieldWriter.dump( "taun", taun ) ;

	//	fieldWriter.dump( "lambda", offline.grains().stresses ) ;
}

int main( int argc, const char* argv[] ) {

	const char * base_dir = "out" ;
	unsigned frame = 0 ;

	bool all = false ;

	bool particles = false ;

	for( int i = 1 ; i < argc ; ++i )
	{
		if( argv[i][0] == '-' ){
			switch(argv[i][1]) {
			case 'a':
				all = true ;
				break ;
			case 'p':
				particles = true ;
				break ;
			case 'n':
				if( ++i == argc ) break ;
				frame = d6::to_uint( argv[i] ) ;
				break ;
			}
		} else {
			base_dir = argv[i] ;
		}
	}

	d6::Offline offline( base_dir ) ;

	unsigned cur_frame = frame ;

	do {
		if(! offline.load_frame( cur_frame ) )
			return all?0:1 ;

		dump_frame( offline, particles, base_dir, cur_frame++ ) ;

	} while( all ) ;

	return (frame == cur_frame) ? 1 : 0 ;

}
