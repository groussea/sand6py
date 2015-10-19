
#include "utils/Config.hh"

#include <cstring>
#include <iostream>

int main( int argc, const char* argv[] )
{
	d6::Config config ;

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
						std::cerr << "Error reading file " << argv[i] << std::endl ;
					}
					break ;
				}
			}
		}
	}

	config.internalize();

	return 0 ;
}
