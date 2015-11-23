#include "utils/Config.hh"
#include "simu/Simu.hh"

#include "utils/Log.hh"
#include "utils/File.hh"

#include <cstring>

int main( int argc, const char* argv[] )
{
	d6::Config config ;

	const char * base_dir = "out" ;

	for( int i = 1 ; i < argc ; ++i )
	{
		if( argv[i][0] == '-' ){
			if( std::strlen(argv[i]) > 2 ) {
				if( ++i == argc ) break ;
				config.from_string( argv[i-1]+1, argv[i] ) ;
			} else {
				switch(argv[i][1]) {
				case 'i':
					if( ++i == argc ) break ;
					if( !config.from_file(argv[i]) ) {
						d6::Log::Error() << "Error reading file " << argv[i] << std::endl ;
					}
					break ;
				case 'v':
					if( ++i == argc ) break ;
					d6::Log::Config::get().setLevel( argv[i] ) ;
				}
			}
		} else {
			base_dir = argv[i] ;
		}
	}

	config.dump( d6::FileInfo(base_dir).filePath("config") );
	config.internalize();

	d6::Simu( config, base_dir ).run() ;

	return 0 ;
}
