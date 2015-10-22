#ifndef D6_OFFLINE_HH
#define D6_OFFLINE_HH

#include "geo/geo.fwd.hh"
#include "geo/Particles.hh"

namespace d6 {

struct Phase ;

class Offline
{
public:

	explicit Offline( const char* base_dir ) ;
	~Offline() ;

	bool load_frame( unsigned frame_nb ) ;

	const Particles &particles() const { return m_particles ; }

	const MeshType &mesh() const { return *m_mesh ; }
	const Phase &phase() const { return *m_phase ; }

private:
	const char* m_base_dir ;

	Particles  m_particles ;
	MeshType*  m_mesh ;
	Phase*     m_phase ;
} ;

} //d6

#endif
