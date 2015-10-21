#ifndef D6_SIMU_HH
#define D6_SIMU_HH

#include "geo/geo.fwd.hh"
#include "simu/DynParticles.hh"

namespace d6 {

class Config ;
struct Phase ;

class Simu {

public:
	typedef Grid MeshImpl ;
	typedef MeshBase< Grid > MeshType ;

	explicit Simu( const Config& config ) ;
	~Simu() ;

	void run() ;
	void step() ;
	void dump(unsigned frame) const ;

	MeshType& mesh() { return * m_mesh ;  }
	const MeshType& mesh() const { return * m_mesh ;  }

private:

	class Fields ;

	Simu( const Simu& ) ;
	Simu& operator=( const Simu& ) ;

	const Config& m_config ;

	DynParticles  m_particles ;
	MeshType*  m_mesh ;
	Phase*     m_phase ;
};

} //d6

#endif
