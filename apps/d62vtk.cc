
#include "utils/string.hh"

#include "visu/Offline.hh"
#include "visu/VTKParticlesWriter.hh"
#include "visu/VTKFieldWriter.hh"

#include "simu/Phase.hh"


int main( int argc, const char* argv[] ) {

	const char * base_dir = "out" ;
	unsigned frame = 0 ;

	for( int i = 1 ; i < argc ; ++i )
	{
		if( argv[i][0] == '-' ){
			switch(argv[i][1]) {
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

	if(! offline.load_frame( frame ) )
		return 1 ;

	d6::VTKParticlesWriter particlesWriter( base_dir, offline.particles() ) ;
	particlesWriter.startFile( "particles", frame ) ;
	particlesWriter.dump_all() ;

	d6::VTKFieldWriter fieldWriter( base_dir, offline.mesh() ) ;
//	fieldWriter.setMode( d6::VTKWriter::Ascii );
	fieldWriter.startFile( "fields", frame ) ;
	fieldWriter.dump(    "phi", offline.grains().fraction ) ;
	fieldWriter.dump(      "u", offline.grains().velocity ) ;
	fieldWriter.dump( "lambda", offline.grains().stresses ) ;

	return 0 ;

}
