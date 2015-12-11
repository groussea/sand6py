#include "utils/Config.hh"
#include "simu/Simu.hh"

#include "utils/Log.hh"
#include "utils/File.hh"

#include <cstring>

#include <fenv.h>

namespace d6 {
	extern const char* g_git_branch ;
	extern const char* g_git_commit ;
	extern const char* g_timestamp  ;
}

int main( int argc, const char* argv[] )
{
	feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);

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

	std::string info = d6::arg( d6::arg3("%1 %3 on %2 [%4]", argv[0], d6::g_git_branch, d6::g_git_commit ), d6::g_timestamp ) ;
	d6::Log::Info() << "This is " << info << std::endl ;

	d6::FileInfo outDir ( base_dir ) ;
	if( !outDir.exists() ) outDir.makeDir() ;
	config.dump( outDir.filePath("config"), info.c_str() );
	config.internalize();

	d6::Log::Verbose() << "Re = " << config.viscosity << std::endl ;

	d6::Simu( config, base_dir ).run() ;

	return 0 ;
}
