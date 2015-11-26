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

	virtual void init( const Params& ) {
		center = Vec( .75*m_config->box[0], .75*m_config->box[1], .5*m_config->box[2] ) ;
		radius = .125 * m_config->box[0] ;
	}

	Scalar particle_density( const Vec &x ) const {
		return ( x[0] > center[0] && x[1] > center[1]
				&& (x - center).squaredNorm() > radius * radius ) ? 1. : 0. ;
	}
};

struct TowerScenar : public Scenario {
	Vec center ;
	Scalar radius ;

	Scalar volMass ;

	Scalar h0   ;
	Scalar hvel ;
	Scalar zvel ;
	Scalar avel ;

	Scalar particle_density( const Vec &x ) const {
		return ( std::fabs( x[0] - center[0] ) < radius
//		        && std::fabs( x[1] - center[1] ) < radius
				)
			   ? 1. : 0. ;
	}

	virtual void init( const Params& params ) {
		center = Vec( .5*m_config->box[0], .5*m_config->box[1], .5*m_config->box[2] ) ;
		radius = .125 * m_config->box[0] ;

		volMass = scalar_param( params,   "vm", Units::VolumicMass, 1.5* m_config->units().R ) ;
		h0      = scalar_param( params,   "h0", Units::Length     , 1. ) ;
		hvel    = scalar_param( params, "hvel", Units::Velocity   , 1. ) ;
		avel    = scalar_param( params, "avel", Units::Frequency  , 0. ) ;

		const Scalar t = h0 / hvel ;
		zvel = .5 * m_config->gravity.norm() * t ;
	}

	void add_rigid_bodies( std::vector< RigidBody >& rbs ) const
	{
		LevelSet::Ptr ls = LevelSet::make_sphere() ;
		ls->scale(.5*.125*m_config->box[0]).set_origin( center - h0*Vec(1,1,0)/M_SQRT2 ) ;

		rbs.emplace_back( ls, volMass );
		rbs.back().set_velocity( Vec(hvel/M_SQRT2, hvel/M_SQRT2, zvel), Vec(avel,0,0) ) ;
	}

	void update( Simu& simu, Scalar /*time*/ ) const
	{
		for( RigidBody& rb: simu.rigidBodies() ) {
			rb.integrate_gravity( m_config->dt(), m_config->gravity );
		}

	}
};

struct RbPlaneTestScenar : public Scenario {
	Scalar particle_density( const Vec &x ) const {
		return ( x[2] >  .5*m_config->box[2] ) ? 1. : 0. ;
	}

	void add_rigid_bodies( std::vector< RigidBody >& rbs ) const
	{
		LevelSet::Ptr ls = LevelSet::make_plane() ;
		ls->set_origin( .5 * m_config->box - Vec(0,0,.25*m_config->box[2]) ) ;
		ls->set_rotation( Vec(1,0,0), M_PI/8 ) ;

		rbs.emplace_back( ls, 1. );
		rbs.back().set_velocity( Vec(0,0,1.e-1), Vec(0,0,0) ) ;
	}
};

struct ImpactScenar : public Scenario {

	Scalar particle_density( const Vec &x ) const {
		return ( x[2] <  .5*m_config->box[2] ) ? 1. : 0. ;
	}

	virtual void init( const Params& params ) {
		volMass = scalar_param( params, "vm", Units::VolumicMass, 1.5*m_config->units().R ) ;
		zvel = scalar_param( params, "zvel", Units::Velocity, 0. ) ;
		avel = scalar_param( params, "avel", Units::Frequency, 0. ) ;
	}

	void add_rigid_bodies( std::vector< RigidBody >& rbs ) const
	{
		LevelSet::Ptr ls = LevelSet::make_sphere() ;
		ls->scale(.125*m_config->box[0]).set_origin( .5 * m_config->box + Vec(0,0,.25*m_config->box[2]) ) ;

		rbs.emplace_back( ls, volMass );
		rbs.back().set_velocity( Vec(0,0,-zvel), Vec(avel,0,0) ) ;
	}

	void update( Simu& simu, Scalar /*time*/ ) const
	{
		for( RigidBody& rb: simu.rigidBodies() ) {
			rb.integrate_gravity( m_config->dt(), m_config->gravity );
		}

	}

private:
	Scalar volMass ;
	Scalar zvel ;
	Scalar avel ;
};



// Factories & stuff

struct DefaultScenarioFactory : public ScenarioFactory
{
	std::unique_ptr< Scenario > make( const std::string & str ) const
	{
		if( str == "rayleigh")
			return std::unique_ptr< Scenario >( new RayleighScenar() ) ;
		if( str == "collapse")
			return std::unique_ptr< Scenario >( new CollapseScenar() ) ;
		if( str == "bridson")
			return std::unique_ptr< Scenario >( new BridsonScenar() ) ;
		if( str == "tower")
			return std::unique_ptr< Scenario >( new TowerScenar() ) ;
		if( str == "rb_plane_test")
			return std::unique_ptr< Scenario >( new RbPlaneTestScenar() ) ;
		if( str == "impact")
			return std::unique_ptr< Scenario >( new ImpactScenar() ) ;

		return std::unique_ptr< Scenario >( new BedScenar() ) ;
	}
} ;

struct ScenarioBuilder {

	static ScenarioBuilder& instance() {
		static ScenarioBuilder s_instance ;
		return s_instance ;
	}

	void add( const ScenarioFactory& factory ) {
		m_factories.push_back( &factory );
	}

	std::unique_ptr< Scenario > make( const std::string& str ) const {
		std::unique_ptr< Scenario > ptr ;
		for( unsigned k = m_factories.size(); k ; --k ) {
			ptr = m_factories[k-1]->make( str ) ;
			if( ptr )
				break ;
		}
		assert( ptr ) ;
		return ptr ;
	}

private:
	ScenarioBuilder()
	{
		add(m_defaultFactory) ;
	}

	DefaultScenarioFactory m_defaultFactory ;
	std::vector< const ScenarioFactory* > m_factories ;
};

std::unique_ptr< Scenario > Scenario::parse( const Config& config )
{
	std::istringstream in( config.scenario ) ;
	std::string line ;
	std::vector< std::string > tok ;

	in >> line ;

	std::unique_ptr< Scenario > pScenar = ScenarioBuilder::instance().make( canonicalize(line) ) ;
	pScenar->m_config = &config ;

	Params params ;
	while( in >> line ) {
		tok.clear() ;
		split( line, ":", tok );
		if( tok.size() == 2 ) {
			params[canonicalize(tok[0])] = canonicalize(tok[1]) ;
		}
	}

	pScenar->init( params ) ;

	return pScenar ;
}

void Scenario::register_factory( const ScenarioFactory& factory )
{
	ScenarioBuilder::instance().add( factory ) ;
}

Scalar Scenario::scalar_param(const Params& params, const std::string& key, Units::Unit unit, Scalar def ) const
{
	Scalar s = def ;
	Params::const_iterator it = params.find(key) ;
	if( it != params.end() ) {
		cast( it->second, s ) ;
	}

	return s * m_config->units().fromSI( unit ) ;

}

} //d6
