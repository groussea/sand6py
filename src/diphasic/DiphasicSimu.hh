#ifndef D6_DIPHASIC_SIMU_HH
#define D6_DIPHASIC_SIMU_HH

#include "simu/Simu.hh"

#include "mono/PhaseMeshes.hh"

namespace d6 {

class DiphasicSimu : public Simu {

public:
	DiphasicSimu( const Config& config, const char* base_dir ) ;
	virtual ~DiphasicSimu() ;

	const PhaseMeshes& meshes() const { return m_meshes ;  }

private:
	PhaseMeshes  m_meshes ;
	std::unique_ptr<Phase>     m_grains ;

};


} //d6


#endif
