#include "Scenario.hh"

#include "Simu.hh"
#include "RigidBody.hh"

#include "geo/LevelSet.hh"

#include "utils/Config.hh"
#include "utils/string.hh"
#include "utils/Log.hh"

namespace d6 {

// Default scenars

struct RayleighScenar : public Scenario {
	Scalar particle_density( const Vec &x ) const override {
		return ( x[1] >  .5*m_config->box[1] ) ? 1. : 0. ;
	}
};

struct BedScenar : public Scenario {
	Scalar particle_density( const Vec &x ) const override {
		return ( x[1] <  .5*m_config->box[1] ) ? 1. : 0. ;
	}
};

struct CollapseScenar : public Scenario {
	Scalar particle_density( const Vec &x ) const override {
		return ( x[0] > .75*m_config->box[0] ) ? 1. : 0. ;
	}

		virtual void init( const Params& params ) override {
	  l0 = scalar_param( params,   "l0", Units::None, .25 ) ;
		}

		private:
		Scalar l0 ;
};

// Factories & stuff

std::unique_ptr< Scenario > DefaultScenarioFactory::make( const std::string & str ) const
{
	if( str == "rayleigh")
		return std::unique_ptr< Scenario >( new RayleighScenar() ) ;
	if( str == "collapse")
		return std::unique_ptr< Scenario >( new CollapseScenar() ) ;

	return std::unique_ptr< Scenario >( new BedScenar() ) ;
}


} //d6
