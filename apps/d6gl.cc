#include "GLViewer.hh"

#include "utils/string.hh"
#include "visu/Offline.hh"

#include <QApplication>

int main( int argc, char* argv[] ) {

	const char * base_dir = "out" ;
	unsigned frame = 0 ;

	bool velCut = false ;
	bool discs  = false ;
	float grainSizeFactor = 1 ;

	QApplication application(argc,argv);

	unsigned nSamples = 0 ;

	for( int i = 1 ; i < argc ; ++i )
	{
		if( argv[i][0] == '-' ){
			switch(argv[i][1]) {
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
			}
		} else {
			base_dir = argv[i] ;
		}
	}

	d6::Offline offline( base_dir ) ;

	QGLFormat fmt = QGLFormat::defaultFormat() ;
	if( nSamples > 0 ) {
		fmt.setSampleBuffers( true );
		fmt.setSamples( 8 );
	}
	d6::GLViewer viewer( fmt, offline, nSamples );

	if( velCut ) {
		viewer.grainsRenderer().cutAndColorVelocities() ;
	} else if ( discs ) {
		viewer.grainsRenderer().useDiscs() ;
	}
	viewer.grainsRenderer().setGrainSizeFactor( grainSizeFactor );

	viewer.set_frame( frame );

	viewer.setWindowTitle("D6 glViewer");

	// Make the viewer window visible on screen.
	viewer.show();

	// Run main loop.
	return application.exec();
}

