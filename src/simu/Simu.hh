#ifndef D6_SIMU_HH
#define D6_SIMU_HH

#include "geo/geo.fwd.hh"
#include "simu/DynParticles.hh"
#include "simu/PhaseSolver.hh"

namespace d6 {

struct Config ;
struct Phase ;
class RigidBody ;

class Simu {

public:
	typedef Grid MeshImpl ;
	typedef MeshBase< Grid > MeshType ;

	explicit Simu( const Config& config, const char* base_dir ) ;
	~Simu() ;

	void run() ;
	void step() ;
	
	void dump_fields(unsigned frame) const ;
	void dump_particles(unsigned frame) const ;

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

private:

	class Fields ;

	Simu( const Simu& ) ;
	Simu& operator=( const Simu& ) ;

	const Config& m_config ;
	const char* m_base_dir ;

	DynParticles  m_particles ;

	std::vector< RigidBody   > m_rigidBodies ;

	MeshType*  m_mesh ;
	Phase*     m_grains ;
	std::vector< TensorField > m_rbStresses  ;

	PhaseSolver m_solver ;
};

} //d6

#endif
