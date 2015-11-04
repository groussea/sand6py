#ifndef D6_SCENARIO_HH
#define D6_SCENARIO_HH

#include "utils/scalar.hh"
#include "geo/Expr.hh"

#include <memory>
#include <string>
#include <unordered_map>

namespace d6 {

class Config ;

class Scenario {

public:

	static std::unique_ptr< Scenario > make( const Config& config ) ;

	struct ParticleGenerator : public Expr<Scalar> {
		const Scenario& scenario ;
		ParticleGenerator( const Scenario& s ) : scenario(s) {}
		Scalar operator() ( const Vec&  x ) const { return scenario.particle_density(x) ; }
	};
	ParticleGenerator generator() const { return ParticleGenerator(*this) ; }

	virtual Scalar particle_density( const Vec &x ) const = 0 ;

protected:
	typedef std::unordered_map< std::string, std::string > Params ;

	Scenario() : m_config(0) {}

	virtual void set_params( const Params& /*params*/ ) {}

	const Config* m_config ;

private:
	static std::unique_ptr< Scenario > make( const std::string& str ) ;
};

} //d6

#endif
