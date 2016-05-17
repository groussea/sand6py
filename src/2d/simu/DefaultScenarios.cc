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
		return ( x[1] >  .5*box()[1] &&
				(x - .5*box()).squaredNorm() > std::pow(box()[0]/64,2) ) ? 1. : 0. ;
	}
};

struct Sedimentation : public Scenario {
	Scalar particle_density( const Vec &x ) const override {
		return 	( x[1] <  .9*box()[1] ) ? .3 : 0 ;
	}
};

struct FallingBallScenar : public Scenario {
	Scalar particle_density( const Vec &x ) const override {
		return (x - (.5*box() + Vec(0,.375*box()[1]) )).squaredNorm() < std::pow(box()[0]/4,2)  ? 1. : 0. ;
	}
};

struct BedScenar : public Scenario {
	Scalar particle_density( const Vec &x ) const override {
		return ( x[1] <  .5*m_config->box[1] ) ? 1. : 0. ;
	}
};

struct CollapseScenar : public Scenario {
	Scalar particle_density( const Vec &x ) const override {
		return ( x[0] > (1-l0)*m_config->box[0] &&  x[1] < .75*m_config->box[1] ) ? 1. : 0. ;
	}

	virtual void init( const Params& params ) override {
		l0 = scalar_param( params,   "l0", Units::None, .25 ) ;
	}

private:
	Scalar l0 ;
};

struct HeapScenar : public Scenario {
	Scalar particle_density( const Vec &x ) const override {
		return ( x[0] > (.375-.175/2)*box()[0] &&  x[0] < (.375+.175/2)*box()[0]
				&& x[1] < .75*box()[1]) ? 1. : 0. ;
	}
};

struct PlaneTestScenar : public Scenario {
	Scalar particle_density( const Vec &x ) const override {
		return ( x[1] >  .5*box(1) ) ? 1. : 0. ;
	}

	void add_rigid_bodies( std::vector< RigidBody >& rbs ) const override
	{
		LevelSet::Ptr ls = LevelSet::make_plane() ;
		ls->set_origin( .5 * m_config->box - Vec(0,.25*m_config->box[1]) ) ;
		ls->set_rotation( M_PI/8 ) ;

		rbs.emplace_back( ls, 1. );
		rbs.back().set_velocity( Vec(0,1.e-1), 0 ) ;
	}
};

struct ImpactScenar : public Scenario {

	Scalar particle_density( const Vec &x ) const override {
		return ( x[1] <  1./3.*m_config->box[1] ) ? 1. : 0. ;
	}

	virtual void init( const Params& params ) {
		volMass = scalar_param( params, "vm", Units::VolumicMass, 1.5*m_config->units().R ) ;
		zvel = scalar_param( params, "zvel", Units::Velocity, 0. ) ;
		avel = scalar_param( params, "avel", Units::Frequency, 0. ) ;
		d = scalar_param( params, "d", Units::None, 0.25 ) ;
	}

	void add_rigid_bodies( std::vector< RigidBody >& rbs ) const override
	{
		LevelSet::Ptr ls = LevelSet::make_sphere() ;
		ls->scale( radius() ).set_origin( .5 * m_config->box + Vec(0,.25*m_config->box[1]) ) ;

		rbs.emplace_back( ls, volMass );
		rbs.back().set_velocity( Vec(0,-zvel), avel ) ;
	}

	Scalar radius() const {
		return d/2*m_config->box[0] ;
	}

	void update( Simu& simu, Scalar /*time*/, Scalar dt ) const override
	{
		for( RigidBody& rb: simu.rigidBodies() ) {
			rb.integrate_gravity( dt, m_config->gravity );
		}

	}

private:
	Scalar volMass ;
	Scalar zvel ;
	Scalar avel ;
	Scalar d ;
};

struct SiloScenar : public Scenario {
	Scalar particle_density( const Vec &x ) const override {
		return ( x[1]-1 >  .5*m_config->box[1] ) ? 1. : 0. ;
	}

	void add_rigid_bodies( std::vector< RigidBody >& rbs ) const override
	{
		const Scalar a = 0.15 ;
		const Scalar L = m_config->box[0] * (1 - a) / 2 ;

		LevelSet::Ptr ls = LevelSet::make_cylinder( L ) ;
		ls->set_origin( Vec( L/2,.5*m_config->box[1]) ) ;
		ls->set_rotation( M_PI/2 ) ;

		LevelSet::Ptr ls2 = LevelSet::make_cylinder( L ) ;
		ls2->set_origin( Vec( m_config->box[0]-L/2,.5*m_config->box[1]) ) ;
		ls2->set_rotation( M_PI/2 ) ;

		rbs.emplace_back( ls2, 1. );
		rbs.emplace_back( ls, 1. );
	}
};

struct TowerScenar : public Scenario {

	Scalar particle_density( const Vec &x ) const override {
		return ( x[0] >= .25*box(0) && x[0] <= .375*box(0) ) ? 1. : 0. ;
	}

	virtual void init( const Params& params ) {
		volMass = scalar_param( params, "vm", Units::VolumicMass, 1.5*m_config->units().R ) ;
		hvel = scalar_param( params, "hvel", Units::Velocity, 1 ) ;
		avel = scalar_param( params, "avel", Units::Frequency, 0. ) ;
		d = scalar_param( params, "d", Units::None, 0.25 ) ;
	}

	void add_rigid_bodies( std::vector< RigidBody >& rbs ) const override
	{
		LevelSet::Ptr ls = LevelSet::make_sphere() ;
		ls->scale( radius() ).set_origin( Vec(0,.375*box(1))  ) ;

		const Scalar t = .25*box(0) / hvel ;
		const Scalar zvel = .5 * m_config->gravity.norm() * t ;

		rbs.emplace_back( ls, volMass );
		rbs.back().set_velocity( Vec(hvel, zvel), avel ) ;
	}

	Scalar radius() const {
		return d/2*m_config->box[1] ;
	}

	void update( Simu& simu, Scalar /*time*/, Scalar dt ) const override
	{
		for( RigidBody& rb: simu.rigidBodies() ) {
			rb.integrate_gravity( dt, m_config->gravity );
		}

	}

private:
	Scalar volMass ;
	Scalar hvel ;
	Scalar avel ;
	Scalar d ;
};

// Factories & stuff

std::unique_ptr< Scenario > DefaultScenarioFactory::make( const std::string & str ) const
{
	if( str == "rayleigh")
		return std::unique_ptr< Scenario >( new RayleighScenar() ) ;
	if( str == "ball")
		return std::unique_ptr< Scenario >( new FallingBallScenar() ) ;
	if( str == "collapse")
		return std::unique_ptr< Scenario >( new CollapseScenar() ) ;
	if( str == "planetest")
		return std::unique_ptr< Scenario >( new PlaneTestScenar() ) ;
	if( str == "impact")
		return std::unique_ptr< Scenario >( new ImpactScenar() ) ;
	if( str == "silo")
		return std::unique_ptr< Scenario >( new SiloScenar() ) ;
	if( str == "tower")
		return std::unique_ptr< Scenario >( new TowerScenar() ) ;
	if( str == "sedim")
		return std::unique_ptr< Scenario >( new Sedimentation() ) ;
	if( str == "heap")
		return std::unique_ptr< Scenario >( new HeapScenar() ) ;

	return std::unique_ptr< Scenario >( new BedScenar() ) ;
}


} //d6
