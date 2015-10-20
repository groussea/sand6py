#include "Config.hh"

#include "string.hh"
#include "Log.hh"

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
	box(1,1,1), res(10,10,10),
	volMass( 1.5e3 ),
	viscosity( 1.e-3 ),
	gravity( 0, 0, -9.81 )
{
}

void Config::internalize()
{
	m_units.setTypical( box.minCoeff(), gravity.norm(), volMass );

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

} //ns hyb2d

