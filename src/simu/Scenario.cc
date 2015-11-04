#include "Scenario.hh"

#include "Config.hh"

#include "utils/string.hh"
#include "utils/Log.hh"

namespace d6 {

struct RayleighScenar : public Scenario {
	Scalar particle_density( const Vec &x ) const {
		return ( x[2] >  .5*m_config->box[2] ) ? 1. : 0. ;
	}
};

struct BedScenar : public Scenario {
	Scalar particle_density( const Vec &x ) const {
		return ( x[2] <  .5*m_config->box[2] ) ? 1. : 0. ;
	}
};

struct CollapseScenar : public Scenario {
	Scalar particle_density( const Vec &x ) const {
		return ( x[0] > .75*m_config->box[0] ) ? 1. : 0. ;
	}
};
struct BridsonScenar : public Scenario {
	Vec center ;
	Scalar radius ;

	virtual void set_params( const Params& params ) {
		Scenario::set_params( params ) ;
		center = Vec( .75*m_config->box[0], .75*m_config->box[1], .5*m_config->box[2] ) ;
		radius = .125 * m_config->box[0] ;
	}

	Scalar particle_density( const Vec &x ) const {
		return ( x[0] > center[0] && x[1] > center[1]
				&& (x - center).squaredNorm() > radius * radius ) ? 1. : 0. ;
	}
};

std::unique_ptr< Scenario > Scenario::make( const std::string& str )
{
	if( str == "rayleigh")
		return std::unique_ptr< Scenario >( new RayleighScenar() ) ;
	if( str == "collapse")
		return std::unique_ptr< Scenario >( new CollapseScenar() ) ;
	if( str == "bridson")
		return std::unique_ptr< Scenario >( new BridsonScenar() ) ;

	return std::unique_ptr< Scenario >( new BedScenar() ) ;
}

std::unique_ptr< Scenario > Scenario::make( const Config& config )
{

	std::istringstream in( config.scenario ) ;
	std::string line ;
	std::vector< std::string > tok ;

	in >> line ;

	std::unique_ptr< Scenario > pScenar = make(canonicalize(line)) ;
	pScenar->m_config = &config ;

	Params params ;
	while( in >> line ) {
		tok.clear() ;
		split( line, ":", tok );
		if( tok.size() == 2 ) {
			params[canonicalize(tok[0])] = canonicalize(tok[1]) ;
		}
	}

	pScenar->set_params( params ) ;

	return pScenar ;
}

} //d6
