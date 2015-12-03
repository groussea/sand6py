#include "Config.hh"

#include "utils/string.hh"
#include "utils/Log.hh"

#include <cmath>

#include <fstream>
#include <regex>

namespace d6 {

template< typename dest_t >
void rescale( dest_t &, Scalar ) { }

template< >
void rescale( Scalar & src, Scalar s ) { src *= s ; }
template< int D >
void rescale( Eigen::Matrix<Scalar, D, 1> &src, Scalar s ) { src *= s ; }


Config::Config() :
	fps(240), substeps(1), nFrames( 1 ),
	box(1,1,1), res(10,10,10), nSamples(2),
	volMass( 1.5e3 ),
	viscosity( 1.e-3 ),
	gravity( 0, 0, -9.81 ),
	phiMax(1), mu(0),
	delta_mu( 0 ), I0( 0.4 ), grainDiameter( 1.e-3 ),
	muRigid( 0.5 ),
	cohesion(0), cohesion_decay(0),
	anisotropy( 0 ),
	enforceMaxFrac( false ),
	boundary("cuve"),
	output( true )
{
}

void Config::internalize()
{
	m_units.setTypical( box.minCoeff()/res.maxCoeff(), gravity.norm(), volMass );

	#define CONFIG_FIELD( name, type, u ) rescale( name, m_units.fromSI( u ) ) ;
	EXPAND_CONFIG
	#undef CONFIG_FIELD
}

bool Config::from_string(const std::string& key, const std::string &val)
{
	std::istringstream iss ( val ) ;
	return from_string( key, iss ) ;
}

bool Config::from_string(const std::string& key, std::istringstream &val )
{
	bool f = false ;
	#define CONFIG_FIELD( name, type, s ) \
		if(key == D6_stringify(name)) { cast(val, name) ; f = true ; }
	EXPAND_CONFIG
	 #undef CONFIG_FIELD

	if( !f ) Log::Warning() << "Warning: '" << key << "' is not a valid config field" << std::endl;
	return f ;
}

bool Config::from_file(const std::string &file_name)
{
	std::ifstream in( file_name ) ;
	if(!in)
		return false ;

	std::string line, key ;

	while( std::getline( in, line ))
	{
		std::istringstream iss ( line ) ;
		if ( (iss >> key) && key[0] != '#' )
		{
			from_string( key, iss ) ;
		}
	}

	return true ;
}

bool Config::dump(const std::string &file_name, const char *comment ) const
{
	std::ofstream out( file_name ) ;
	if(!out)
		return false ;

	if(comment) {
		out << "# " << comment << std::endl ;
	}

	#define CONFIG_FIELD( name, type, s ) \
		out << D6_stringify(name) << "\t" ; d6::dump( out, name ) ; out << "\n" ;
	EXPAND_CONFIG
	 #undef CONFIG_FIELD

	return true ;
}

} //ns hyb2d

