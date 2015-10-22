
#include "Offline.hh"

#include "simu/Phase.hh"

#include "geo/Grid.hh"
#include "utils/Log.hh"
#include "utils/File.hh"

#include "geo/Particles.io.hh"
#include <boost/archive/binary_iarchive.hpp>

#include <bogus/Core/Eigen/EigenSerialization.hpp>

namespace d6 {

Offline::Offline(const char *base_dir)
	: m_base_dir( base_dir )
{
	m_mesh = new MeshImpl( Vec::Ones(), Vec3i::Ones() ) ;
	m_phase = new Phase( *m_mesh ) ;
}

Offline::~Offline(){
	delete m_mesh ;
	delete m_phase ;
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
		ia >> *m_phase ;
	}
	// Particles
	{
		std::ifstream ifs( dir.filePath("particles") );
		boost::archive::binary_iarchive ia(ifs);
		ia >> m_particles ;
	}

	return true ;

}


} //d6
