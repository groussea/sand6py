#include "GLViewer.hh"

#include "utils/string.hh"
#include "visu/Offline.hh"

#include <QApplication>

static void usage( const char *name )
{
	std::cout << "Usage: " << name
			  << " [sim_dir=out] [options] "
			  << "\nOpenGL viewer, based on libQGLViewer."
			  << "\n Press 'H' for help."
			  << "\n\n" ;

	std::cout << "Options:\n"
			  << "-? \t Display this help message and exit\n"
			  << "-n frame_id \t Jump to frame frame_id \n"
			  << "-s n \t Render n samples per particle \n"
			  << "-v   \t Orthographic velocity view \n"
			  << "-d   \t Render discs instead of grains \n"
			  << "-g size \t Multiplicative factor for the size of rendered samples\n"
			  << "-w width \t Specify viewport width\n"
			  << "-h height \t Specify viewport height\n"
			  << "-l dir \t Specify light direction \n"
			  << std::endl ;
}

int main( int argc, char* argv[] ) {

	const char * base_dir = "out" ;
	unsigned frame = 0 ;

	bool velCut = false ;
	bool discs  = false ;
	float grainSizeFactor = 1 ;

	unsigned width  = 1280 ;
	unsigned height = 720 ;
	Eigen::Vector3f lightDir (-.5,-1,1.) ;

	QApplication application(argc,argv);

	unsigned nSamples = 0 ;

	for( int i = 1 ; i < argc ; ++i )
	{
		if( argv[i][0] == '-' ){
			switch(argv[i][1]) {
			case '?':
				usage( argv[0]) ;
				return 0;
			case 'n':
				if( ++i == argc ) break ;
				frame = d6::to_uint( argv[i] ) ;
				break ;
			case 's':
				if( ++i == argc ) break ;
				nSamples = d6::to_uint( argv[i] ) ;
				break;
			case 'v':
				velCut = true ;
				break;
			case 'd':
				discs = true ;
				break;
			case 'g':
				if( ++i == argc ) break ;
				grainSizeFactor = d6::to_float( argv[i] ) ;
				break;
			case 'w':
				if( ++i == argc ) break ;
				width = d6::to_uint( argv[i] ) ;
				break;
			case 'h':
				if( ++i == argc ) break ;
				height = d6::to_uint( argv[i] ) ;
				break;
			case 'l':
				if( ++i == argc ) break ;
								d6::cast( argv[i], lightDir ) ;
				break;
			}
		} else {
			base_dir = argv[i] ;
		}
	}

	d6::Offline offline( base_dir ) ;

	QGLFormat fmt = QGLFormat::defaultFormat() ;
	fmt.setSampleBuffers( true );
	fmt.setSamples( 8 );
	d6::GLViewer viewer( fmt, offline, nSamples, width, height );

	if( velCut ) {
		viewer.grainsRenderer().cutAndColorVelocities() ;
	} else if ( discs ) {
		viewer.grainsRenderer().useDiscs() ;
	}
	viewer.grainsRenderer().setGrainSizeFactor( grainSizeFactor );

	viewer.set_frame( frame );
		viewer.setLightDirection( lightDir ) ;

	viewer.setWindowTitle("D6 glViewer");

	// Make the viewer window visible on screen.
	viewer.show();

	// Run main loop.
	return application.exec();
}

