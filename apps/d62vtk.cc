
#include "utils/string.hh"

#include "visu/Offline.hh"
#include "visu/VTKParticlesWriter.hh"


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
	particlesWriter.dump( frame, d6::VTKParticlesWriter::Volumes );
	particlesWriter.dump( frame, d6::VTKParticlesWriter::Velocities );
	particlesWriter.dump( frame, d6::VTKParticlesWriter::Frames );

	return 0 ;

}
