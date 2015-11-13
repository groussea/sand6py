#include "GLViewer.hh"

#include "utils/string.hh"
#include "visu/Offline.hh"

#include <QApplication>

int main( int argc, char* argv[] ) {

	const char * base_dir = "out" ;
	unsigned frame = 0 ;

	QApplication application(argc,argv);

	bool samples = false ;

	for( int i = 1 ; i < argc ; ++i )
	{
		if( argv[i][0] == '-' ){
			switch(argv[i][1]) {
			case 'n':
				if( ++i == argc ) break ;
				frame = d6::to_uint( argv[i] ) ;
				break ;
			case 's':
				samples = true ;
				break;
			}
		} else {
			base_dir = argv[i] ;
		}
	}

	d6::Offline offline( base_dir ) ;

	d6::GLViewer viewer( offline, samples );
	viewer.set_frame( frame );

	viewer.setWindowTitle("D6 glViewer");

	// Make the viewer window visible on screen.
	viewer.show();

	// Run main loop.
	return application.exec();
}

