
#include "Offline.hh"

#include "mono/Phase.hh"

#include "geo/MeshImpl.hh"

#include "utils/Log.hh"
#include "utils/File.hh"

#include "geo/Particles.io.hh"
#include "geo/LevelSet.io.hh"
#include "geo/Meshes.io.hh"

#include "utils/serialization.hh"

#include <boost/archive/binary_iarchive.hpp>

namespace d6 {

Offline::Offline(const char *base_dir)
	: m_base_dir( base_dir ),
	  m_meshes{ std::unique_ptr<PrimalMesh>(new PrimalMesh( Vec::Ones(), VecWi::Ones(), &m_particles )),
				std::unique_ptr<  DualMesh>(new   DualMesh( Vec::Ones(), VecWi::Ones(), &m_particles ))
			   },
	  m_grains( new Phase( m_meshes ) )
{
	m_config.from_file( FileInfo( base_dir ).filePath( "config" ) ) ;
	m_config.internalize() ;
	Log::Info() << arg( "Frame Dt is %1 (IU) = %2 FPS (SI)", frame_dt(), m_config.fps / m_config.units().toSI( Units::Time )) << std::endl ;;
}

Offline::~Offline(){
}

bool Offline::load_frame(unsigned frame )
{
	FileInfo dir( FileInfo(m_base_dir).filePath( arg("frame-%1", frame ) ) ) ;
	if( !dir.exists() ) {
		Log::Error() << "Directory does not exist: " << dir.path() << std::endl ;
		return false ;
	}

	try {

		// Particles
		{
			std::ifstream ifs( dir.filePath("particles") );
			boost::archive::binary_iarchive ia(ifs);
			ia >> m_particles ;
		}
		// Grid
		{
			std::ifstream ifs( dir.filePath("meshes") );
			boost::archive::binary_iarchive ia(ifs);
			ia >> m_meshes ;
		}
		// Velocity, Stress, Phi
		{
			std::ifstream ifs( dir.filePath("fields") );
			boost::archive::binary_iarchive ia(ifs);
			ia >> *m_grains ;
		}
		// Log
		{
			m_events.clear() ;
			std::ifstream ifs( dir.filePath("log") );
			if(ifs) {
				boost::archive::binary_iarchive ia(ifs);
				ia >> m_events ;
			}
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

	} catch (boost::archive::archive_exception &e) {
		Log::Error() << "Error reading frame data: " << e.what() << std::endl ;
		return false ;
	}

	Log::Info() << "Loaded frame " << frame << std::endl ;

	return true ;

}

} //d6
