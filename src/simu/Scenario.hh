#ifndef D6_SCENARIO_HH
#define D6_SCENARIO_HH

#include "utils/scalar.hh"
#include "geo/Expr.hh"

#include <memory>
#include <string>
#include <unordered_map>

namespace d6 {

class Config ;
class ScenarioFactory ;
class RigidBody ;
class Simu ;

class Scenario {

public:

	static std::unique_ptr< Scenario > parse( const Config& config ) ;
	static void register_factory( const ScenarioFactory& factory ) ;

	virtual ~Scenario() {}

	struct ParticleGenerator : public Expr<Scalar> {
		const Scenario& scenario ;
		ParticleGenerator( const Scenario& s ) : scenario(s) {}
		Scalar operator() ( const Vec&  x ) const { return scenario.particle_density(x) ; }
	};
	ParticleGenerator generator() const { return ParticleGenerator(*this) ; }

	virtual Scalar particle_density( const Vec &x ) const = 0 ;

	virtual void add_rigid_bodies( std::vector< RigidBody >& /*rbs*/ ) const {}
	virtual void update( Simu& /*simu*/, Scalar /*time*/ ) const {}

protected:
	typedef std::unordered_map< std::string, std::string > Params ;

	Scenario() : m_config(0) {}
	virtual void init( const Params& /*params*/ ) {}

	const Config* m_config ;

};

struct ScenarioFactory {
	virtual std::unique_ptr< Scenario > make( const std::string& str ) const = 0 ;

	ScenarioFactory() {}
	ScenarioFactory( const ScenarioFactory& ) = delete ;
	ScenarioFactory& operator= ( const ScenarioFactory& ) = delete ;
};


} //d6

#endif
