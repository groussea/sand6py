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
		return ( x[2] >  .5*m_config->box[2] ) ? 1. : 0. ;
	}
};

struct BedScenar : public Scenario {
	Scalar particle_density( const Vec &x ) const override {
		return ( x[2] <  .5*m_config->box[2] ) ? 1. : 0. ;
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
struct BridsonScenar : public Scenario {
	Vec center ;
	Scalar radius ;

	virtual void init( const Params& ) override {
		center = Vec( .75*m_config->box[0], .75*m_config->box[1], .5*m_config->box[2] ) ;
		radius = .125 * m_config->box[0] ;
	}

	Scalar particle_density( const Vec &x ) const override {
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

	Scalar particle_density( const Vec &x ) const override {
		return (
					( std::fabs( x[0] - center[0] ) < radius/2
				&& std::fabs( x[1] - center[1] ) < radius )
				||	( std::fabs( x[0] - center[0] ) < radius
				&& std::fabs( x[1] - center[1] ) < radius/2 )

				//std::pow( x[1] - center[1], 2 ) + std::pow( x[0] - center[0], 2 ) < radius * radius
				|| std::fabs( x[2] ) < .25* radius
				)
			   ? 1. : 0. ;
	}

	virtual void init( const Params& params ) override {
		center = Vec( .5*m_config->box[0], .5*m_config->box[1], .5*m_config->box[2] ) ;
		radius = .125 * m_config->box[0] ;

		volMass = scalar_param( params,   "vm", Units::VolumicMass, 1.5* m_config->units().R ) ;
		h0      = scalar_param( params,   "h0", Units::Length     , 1. ) ;
		hvel    = scalar_param( params, "hvel", Units::Velocity   , 1. ) ;
		avel    = scalar_param( params, "avel", Units::Frequency  , 0. ) ;

		const Scalar t = h0 / hvel ;
		zvel = .5 * m_config->gravity.norm() * t ;
	}

	void add_rigid_bodies( std::vector< RigidBody >& rbs ) const override
	{
		LevelSet::Ptr ls = LevelSet::make_sphere() ;
		Vec pos = center - h0*Vec(1,1,0)/M_SQRT2 ;
		pos[2] -= .125 * m_config->box[2];
		ls->scale(.5*.125*m_config->box[0]).set_origin( pos ) ;

		rbs.emplace_back( ls, volMass );
		rbs.back().set_velocity( Vec(hvel/M_SQRT2, hvel/M_SQRT2, zvel), Vec(avel,0,0) ) ;
	}

	void update( Simu& simu, Scalar /*time*/, Scalar dt ) const override
	{
		for( RigidBody& rb: simu.rigidBodies() ) {
			rb.integrate_gravity( dt, m_config->gravity );
		}

	}
};

struct RbPlaneTestScenar : public Scenario {
	Scalar particle_density( const Vec &x ) const override {
		return ( x[2] >  .5*m_config->box[2] ) ? 1. : 0. ;
	}

	void add_rigid_bodies( std::vector< RigidBody >& rbs ) const override
	{
		LevelSet::Ptr ls = LevelSet::make_plane() ;
		ls->set_origin( .5 * m_config->box - Vec(0,0,.25*m_config->box[2]) ) ;
		ls->set_rotation( Vec(1,0,0), M_PI/8 ) ;

		rbs.emplace_back( ls, 1. );
		rbs.back().set_velocity( Vec(0,0,1.e-1), Vec(0,0,0) ) ;
	}
};

struct ImpactScenar : public Scenario {

	Scalar particle_density( const Vec &x ) const override {
		return ( x[2] <  1./3.*m_config->box[2] ) ? 1. : 0. ;
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
		ls->scale( radius() ).set_origin( .5 * m_config->box + Vec(0,0,.25*m_config->box[2]) ) ;

		rbs.emplace_back( ls, volMass );
		rbs.back().set_velocity( Vec(0,0,-zvel), Vec(avel,0,0) ) ;
	}

	void update( Simu& simu, Scalar /*time*/, Scalar dt ) const override
	{
		for( RigidBody& rb: simu.rigidBodies() ) {
			rb.integrate_gravity( dt, m_config->gravity );
			Vec vel = rb.velocity() ;
			vel[2] = std::max( vel[2], (radius() - rb.levelSet().origin()[2])/dt ) ;
			rb.set_velocity( vel, rb.angularVelocity() );
		}

	}

	Scalar radius() const {
		return d/2*m_config->box[0] ;
	}

private:
	Scalar volMass ;
	Scalar zvel ;
	Scalar avel ;
	Scalar d ;
};

struct WheelScenar : public Scenario {

	Scalar particle_density( const Vec &x ) const override {
		return ( x[2] <  1./2.*m_config->box[2] ) ? 1. : 0. ;
	}

	virtual void init( const Params& params ) {
		volMass = scalar_param( params, "vm", Units::VolumicMass, 1.5*m_config->units().R ) ;
		xvel = scalar_param( params, "xvel", Units::Velocity, 0. ) ;
		avel = scalar_param( params, "avel", Units::Frequency, 0. ) ;
		torque = scalar_param( params, "torque", Units::Torque, 0. ) ;
		tilt = scalar_param( params, "tilt", Units::None, 0. ) ;
		R = scalar_param( params, "R", Units::None, 0.125 ) ;
		r = scalar_param( params, "r", Units::None, 0.375 ) ;
		angFric = scalar_param( params, "af", Units::None, 0 )
				* m_config->units().fromSI( Units::Torque )
				* m_config->units().fromSI( Units::Time )
				;
	}

	void add_rigid_bodies( std::vector< RigidBody >& rbs ) const override
	{
		LevelSet::Ptr ls = LevelSet::make_torus(r) ;
		ls->scale( radius() )
				.rotate(Vec(1,0,0), (1-tilt) * M_PI_2)
				.set_origin( .5 * m_config->box + Vec(.25*m_config->box[0],-.25*m_config->box[1], radius()*(1+2*r)) ) ;

		rbs.emplace_back( ls, volMass );
		rbs.back().set_velocity( Vec(0,0,-xvel), Vec(0,-avel,0) ) ;
	}

	void update( Simu& simu, Scalar /*time*/, Scalar dt ) const override
	{
		Vec6 forces ;
		forces.head<3>().setZero() ;

		for( RigidBody& rb: simu.rigidBodies() ) {
			forces.tail<3>() = rb.levelSet().rotation()
					* Vec(0,0,torque - angFric * rb.angularVelocity().norm() ) ;
			rb.integrate_gravity( dt, m_config->gravity );
			rb.integrate_forces ( dt, forces );
		}

	}

	Scalar radius() const {
		return R*m_config->box[0] ;
	}

private:
	Scalar volMass ;
	Scalar xvel ;
	Scalar avel ;
	Scalar torque ;
	Scalar tilt ;
	Scalar R ;
	Scalar r ;
	Scalar angFric ;
};

struct SiloScenar : public Scenario {

	constexpr static Scalar s_h = 2 ;

	Scalar particle_density( const Vec &x ) const override {
		return ( x[2] > ( .5*m_config->box[2] + s_h ) //(1-Dbar)*R )
			//	|| x[2] < h*m_config->box[2]
				) ? 1. : 0. ;
	}

	virtual void init( const Params& params ) override {
		W = .5*m_config->box[0] ;
		R = scalar_param( params, "d", Units::None, 0.5 ) * W ;
	}

	void add_rigid_bodies( std::vector< RigidBody >& rbs ) const override
	{
		LevelSet::Ptr ls = LevelSet::make_hole( R/s_h + 1 ) ;
		ls->scale(s_h).set_origin( .5 * m_config->box ) ;
		rbs.emplace_back( ls, 1.e99 );
	}

	void update( Simu& simu, Scalar /*time*/, Scalar /*dt*/ ) const override
	{
		DynParticles &particles = simu.particles() ;
		const Scalar zmin = m_config->box[2] / 3 ;

#pragma omp parallel for
		for( size_t i = 0 ; i < particles.count() ; ++i ) {
			if( particles.geo().centers()(2,i) < zmin ) {
				particles.remove( i );
			}
		}
	}

private:
	Scalar W ;
	Scalar R ;
};

struct HourglassScenar : public Scenario {

	Scalar particle_density( const Vec &x ) const override {
		return ( hg_ls->eval_at( x ) < 0 &&
				 x[2] > .5*m_config->box[2] && x[2] < .8*m_config->box[2]
		) ? 1 : 0 ;
	}

	void init( const Params& params ) override {
		const Scalar D = scalar_param( params, "D", Units::None, 0.9 ) ;
		const Scalar d = scalar_param( params, "d", Units::None, 0.25 ) ;

		const Scalar S = .5*D*m_config->box[1] ;
		const Scalar H = m_config->box[2]/2 / S ;

		hg_ls = LevelSet::make_hourglass( H, d ) ;
		hg_ls->scale( S ).set_origin( .5 * m_config->box ) ;
	}

	void add_rigid_bodies( std::vector< RigidBody >& rbs ) const override
	{
		rbs.emplace_back( hg_ls, 1.e99 );
	}

	mutable LevelSet::Ptr hg_ls ;
};

struct BunnyScenar : public Scenario {

	Scalar particle_density( const Vec &x ) const override {
		return ( bunny_ls->eval_at( x ) < 0 &&
				 x[2] < .5*m_config->box[2]
		) ? 1 : 0 ;
	}

	void init( const Params& params ) override {
		const std::string& meshname = string_param( params, "mesh", "../scenes/bunny.obj") ;
		bunny_ls = LevelSet::from_mesh( meshname.c_str() ) ;

		Eigen::AngleAxis< Scalar > aa( ) ;

		const Scalar S = m_config->box[0] ;

		bunny_ls->scale( S*4 )
				.rotate( Vec(0,1,0), M_PI/4 )
				.rotate( Vec(1,0,0), M_PI/4 )
				.set_origin( S * Vec(.5,.25,-.375) ) ;
		bunny_ls->compute( ) ;
	}

	void add_rigid_bodies( std::vector< RigidBody >& rbs ) const override
	{
		rbs.emplace_back( bunny_ls, 1.e99 );
	}

	void update( Simu& simu, Scalar time, Scalar /*dt*/ ) const override
	{
		const Scalar speed = m_config->units().toSI( Units::Time ) ;

		const Scalar tw = time * speed ;
		Vec vel = Vec::Zero() ;

		if( tw > .25 ) {
			const Scalar t = tw-.25 ;
			if( t < 1 ) {
				vel = Vec(0,0,6*(t - t*t)) * speed * ( .25 * m_config->box[2] );
			}
		}

		for( RigidBody& rb: simu.rigidBodies() ) {
			rb.set_velocity( vel, Vec::Zero() );
		}
	}

	mutable LevelSet::Ptr bunny_ls ;
};

struct WritingScenar : public Scenario {

	Scalar particle_density( const Vec &x ) const override {
		return ( x[2] < .5*m_config->box[2]	) ? 1 : 0 ;
	}

	void init( const Params& params ) override {

		R = scalar_param( params, "d", Units::None, 0.05 ) * m_config->box[0] ;
		H = scalar_param( params, "h", Units::None, 0.5 ) * m_config->box[2] ;
		S = scalar_param( params, "s", Units::None, 1. ) ;

		bezier.resize(2,12) ;
		bezier <<  0.196,  1.046, -0.358, 0.492, 0.609, 0.858, 0.677, 0.737, 1.068, 0.518, 1.630, 0.990,
				   0.140,  0.269,  0.450, 0.561, 0.482, 0.686, 0.630, 0.106, 0.574, -0.25, 0.231, 0.364 ;

		bezier *= m_config->box[0] / 1.4 ;
	}

	void add_rigid_bodies( std::vector< RigidBody >& rbs ) const override
	{
		LevelSet::Ptr cyl_ls ;
		cyl_ls = LevelSet::make_cylinder( H/R ) ;
		cyl_ls->scale(R).set_origin(Vec::Zero()) ;
		rbs.emplace_back( cyl_ls, 1.e99 );
	}

	void update( Simu& simu, Scalar time, Scalar /*dt*/ ) const override
	{
		const unsigned nCurves = bezier.cols() / 4 ;
		const Scalar speed = m_config->units().toSI( Units::Time ) * S ;
		const Scalar trans = .1 ;

		const Scalar tw = time * speed ;
		const unsigned cid = std::floor( tw / ( 1. + trans ) ) ;

		Scalar t = tw - ( 1. + trans )*cid ;
		const bool inbetween = t < trans ;
		t -= trans ;

		Vec vel = Vec::Zero() ;

		if( cid < nCurves && !inbetween ) {
			vel.head<2>() = -3*(1-t)*(1-t)    * bezier.col(0+4*cid)
					+ 3*( (1-t) - 2*t )*(1-t) * bezier.col(1+4*cid)
					+ 3*( 2*(1-t) - t )*t     * bezier.col(2+4*cid)
					+ 3*t*t                   * bezier.col(3+4*cid)
					;
			const Scalar S = 5.e1 ;
			vel[2] = - (H/2-R)*2*S*t*std::exp(-S*t*t) ;
			vel *= speed ;
		}

		for( RigidBody& rb: simu.rigidBodies() ) {
			if( cid < nCurves && inbetween ) {
				Vec origin ;
				origin.head<2>() = bezier.col(4*cid) ;
				origin[2] = .5 * m_config->box[2] + H/2 ;
				rb.move_to( origin );
			}

			rb.set_velocity( vel, Vec::Zero() );
		}

	}

	Eigen::Matrix< Scalar, 2, Eigen::Dynamic > bezier ;
	Scalar R ;
	Scalar H ;
	Scalar S ;

};

struct DiggingScenar : public Scenario {

	Scalar particle_density( const Vec &x ) const override {
		return ( x[2] < .5*m_config->box[2]
		) ? 1 : 0 ;
	}

	void init( const Params& params ) override {
		speedup = scalar_param( params, "s", Units::None, .5 ) ;

		std::string meshname = string_param( params, "mesh", "../scenes/Schmuck2l.obj") ;
		left_hand_ls = LevelSet::from_mesh( meshname.c_str() ) ;
		meshname = string_param( params, "mesh", "../scenes/Schmuck2r.obj") ;
		right_hand_ls = LevelSet::from_mesh( meshname.c_str() ) ;

		const Scalar S = m_config->box[0] ;
		const Scalar H = m_config->box[2] ;
		const Scalar dist = .85 ;

		left_hand_ls->scale( .5*S )
				.rotate( Vec(0,1,0), M_PI_2 )
				.set_origin( Vec(S/2,(1-dist)*S,H) ) ;
		left_hand_ls->compute( ) ;

		right_hand_ls->scale( .5*S )
				.rotate( Vec(0,1,0), -M_PI_2 )
				.rotate( Vec(0,0,1), -M_PI )
				.set_origin( Vec(S/2,dist*S,H) ) ;
		right_hand_ls->compute( ) ;
	}

	void add_rigid_bodies( std::vector< RigidBody >& rbs ) const override
	{
		rbs.emplace_back( left_hand_ls, 1.e99 );
		rbs.emplace_back( right_hand_ls, 1.e99 );
	}

	void update( Simu& simu, Scalar time, Scalar /*dt*/ ) const override
	{
		const Scalar speed = speedup * M_PI_2 * m_config->units().toSI( Units::Time ) ;

		const Scalar tw = time * speed ;
		Vec  vel = Vec::Zero() ;
		Vec avel = Vec::Zero() ;

		if( tw > 1 ) {
			const Scalar t = tw - 1;
			if( t < 1 ) {
				avel = Vec(0,0,0);
				vel = Vec(0, 0, 6*(t-t*t) ) * .5 * m_config->box[2] ;
			}
		}else if( tw > 0 ) {
			const Scalar t = tw;
			if( t < 1 ) {
				avel = Vec(.8*M_PI_2,0,0);
				vel = Vec(0, std::sin(t), -std::cos(t)) * .9 * m_config->box[2] ;
			}
		}

		vel *= speed ;
		avel *= speed ;

		int i = 1 ;
		for( RigidBody& rb: simu.rigidBodies() ) {
			rb.set_velocity( Vec(vel[0], i*vel[1], vel[2]), i * avel );
			i = -i ;
		}
	}

	Scalar speedup ;
	mutable LevelSet::Ptr  left_hand_ls ;
	mutable LevelSet::Ptr right_hand_ls ;
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
		if( str == "silo")
			return std::unique_ptr< Scenario >( new SiloScenar() ) ;
		if( str == "bunny")
			return std::unique_ptr< Scenario >( new BunnyScenar() ) ;
		if( str == "writing")
			return std::unique_ptr< Scenario >( new WritingScenar() ) ;
		if( str == "hourglass")
			return std::unique_ptr< Scenario >( new HourglassScenar() ) ;
		if( str == "wheel")
			return std::unique_ptr< Scenario >( new WheelScenar() ) ;
		if( str == "digging")
			return std::unique_ptr< Scenario >( new DiggingScenar() ) ;

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
std::string Scenario::string_param(const Params &params, const std::string &key, const std::string &def) const
{
	std::string s = def ;
	Params::const_iterator it = params.find(key) ;
	if( it != params.end() ) {
		s = it->second ;
	}
	return s ;

}

} //d6
