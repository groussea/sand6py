#ifndef D6_OFFLINE_HH
#define D6_OFFLINE_HH

#include "geo/geo.fwd.hh"
#include "geo/Particles.hh"

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

private:
	const char* m_base_dir ;

	Particles  m_particles ;
	MeshType*  m_mesh ;
	Phase*     m_grains ;

	std::vector< std::unique_ptr< LevelSet > > m_levelSets ;
} ;

} //d6

#endif
