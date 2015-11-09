
#include "Offline.hh"

#include "simu/Phase.hh"

#include "geo/Grid.hh"
#include "utils/Log.hh"
#include "utils/File.hh"

#include "geo/Particles.io.hh"
#include "geo/LevelSet.io.hh"

#include <boost/archive/binary_iarchive.hpp>

#include <bogus/Core/Eigen/EigenSerialization.hpp>

namespace d6 {

Offline::Offline(const char *base_dir)
	: m_base_dir( base_dir )
{
	m_mesh = new MeshImpl( Vec::Ones(), Vec3i::Ones() ) ;
	m_grains = new Phase( *m_mesh ) ;
}

Offline::~Offline(){
	delete m_mesh ;
	delete m_grains ;
}

bool Offline::load_frame(unsigned frame )
{
	FileInfo dir( FileInfo(m_base_dir).filePath( arg("frame-%1", frame ) ) ) ;
	if( !dir.exists() ) {
		Log::Error() << "Directory does not exist: " << dir.path() << std::endl ;
		return false ;
	}

	// Grid
	{
		std::ifstream ifs( dir.filePath("mesh") );
		boost::archive::binary_iarchive ia(ifs);
		ia >> m_mesh->derived() ;
	}
	// Velocity, Stress, Phi
	{
		std::ifstream ifs( dir.filePath("fields") );
		boost::archive::binary_iarchive ia(ifs);
		ia >> *m_grains ;
	}
	// Particles
	{
		std::ifstream ifs( dir.filePath("particles") );
		boost::archive::binary_iarchive ia(ifs);
		ia >> m_particles ;
	}

	//Objects
	{
		m_levelSets.clear();

		std::ifstream ifs( dir.filePath("objects") );
		boost::archive::binary_iarchive ia(ifs);

		unsigned n = 0 ;
		ia >> n ;
		LevelSet::register_derived(ia) ;
		for( unsigned i = 0 ; i < n ; ++i ) {
			LevelSet* ptr ;
			ia >> ptr ;
			m_levelSets.emplace_back( ptr ) ;
		}

	}

	Log::Info() << "Load frame " << frame << std::endl ;

	return true ;

}


} //d6
