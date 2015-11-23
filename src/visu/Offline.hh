#ifndef D6_OFFLINE_HH
#define D6_OFFLINE_HH


#include "geo/geo.fwd.hh"
#include "geo/Particles.hh"

#include "utils/Config.hh"

#include <memory>

namespace d6 {

struct Phase ;
class LevelSet ;

class Offline
{
public:

	explicit Offline( const char* base_dir ) ;
	~Offline() ;

	bool load_frame( unsigned frame_nb ) ;

	const Particles &particles() const { return m_particles ; }

	const MeshType &mesh() const { return *m_mesh ; }
	const Phase &grains() const { return *m_grains ; }

	const char* base_dir() const
	{
		return m_base_dir ;
	}

	const std::vector< std::unique_ptr< LevelSet > >& levelSets() const
	{
		return m_levelSets ;
	}

	Scalar frame_dt () const {
		return 1./ m_config.fps ;
	}

	const Particles::EventLog& events() const {
		return m_events ;
	}

private:
	const char* m_base_dir ;
	Config m_config ;

	Particles  m_particles ;
	Particles::EventLog  m_events ;

	MeshType*  m_mesh ;
	Phase*     m_grains ;

	std::vector< std::unique_ptr< LevelSet > > m_levelSets ;
} ;

} //d6

#endif
