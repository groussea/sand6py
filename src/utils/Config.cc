#include "Config.hh"

#include "string.hh"

#include <cmath>

#include <iostream>
#include <fstream>
#include <regex>

#define stringify(s) preproc_str(s)
#define preproc_str(s) #s

namespace d6 {

template< typename dest_t >
void rescale( dest_t &, Scalar ) { }

template< >
void rescale( Scalar & src, Scalar s ) { src *= s ; }
template< int D >
void rescale( Eigen::Matrix<Scalar, D, 1> &src, Scalar s ) { src *= s ; }


Config::Config() :
	fps(240), substeps(1),
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

void show_matches(const std::string& in, const std::string& re)
{
	std::smatch m;
	std::regex_search(in, m, std::regex(re));
	if(m.empty()) {
		std::cout << "input=[" << in << "], regex=[" << re << "]: NO MATCH\n";
	} else {
		std::cout << "input=[" << in << "], regex=[" << re << "]: ";
		std::cout << "prefix=[" << m.prefix() << "] ";
		for(std::size_t n = 0; n < m.size(); ++n)
			std::cout << " m[" << n << "]=[" << m[n] << "] ";
		std::cout << "suffix=[" << m.suffix() << "]\n";
	}
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
		if(key == stringify(name)) { cast(val, name) ; f = true ; }
	EXPAND_CONFIG
	 #undef CONFIG_FIELD

	if( !f ) std::cerr << "Warning: '" << key << "' is not a valid config field" << std::endl;
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

