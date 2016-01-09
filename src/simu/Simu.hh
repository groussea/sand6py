#ifndef D6_SIMU_HH
#define D6_SIMU_HH

#include "geo/geo.fwd.hh"

#include "simu/DynParticles.hh"
#include "simu/PhaseSolver.hh"

#include "utils/Stats.hh"

#include <memory>

namespace d6 {

struct Config ;
struct Phase ;
class RigidBody ;
class Scenario ;

class Simu {

public:

	explicit Simu( const Config& config, const char* base_dir ) ;
	~Simu() ;

	//! Runs the simulation
	void run() ;

	MeshType& mesh() { return * m_mesh ;  }
	const MeshType& mesh() const { return * m_mesh ;  }

	std::vector< RigidBody > &rigidBodies ()
	{
		return m_rigidBodies ;
	}
	const std::vector< RigidBody > &rigidBodies () const
	{
		return m_rigidBodies ;
	}

	DynParticles& particles() { return m_particles ; }
	const DynParticles& particles() const { return m_particles ; }

private:
	//! Advances the simulation with time step \p dt
	void step( const Scalar dt ) ;

	// Output
	void dump_fields(unsigned frame) const ;
	void dump_particles(unsigned frame) const ;

	Simu( const Simu& ) = delete ;
	Simu& operator=( const Simu& ) = delete ;


	const Config& m_config ;
	const char* m_base_dir ;

	Stats		  m_stats ;
	DynParticles  m_particles ;
	std::unique_ptr<Scenario>  m_scenario ;

	std::vector< RigidBody   > m_rigidBodies ;

	MeshType*  m_mesh ;
	Phase*     m_grains ;
	// Useful for warm-starting stresses at frictional boundary conditions
	std::vector< TensorField > m_rbStresses  ;

	PhaseSolver m_solver ;
};

} //d6

#endif
