#ifndef D6_DIPHASIC_SIMU_HH
#define D6_DIPHASIC_SIMU_HH

#include "simu/Simu.hh"

#include "DiphasicSolver.hh"

#include "mono/PhaseMeshes.hh"

namespace d6 {

class DiphasicSimu : public Simu {

public:
	DiphasicSimu( const Config& config, const char* base_dir ) ;
	virtual ~DiphasicSimu() ;

	const PhaseMeshes& meshes() const { return m_meshes ;  }

protected:

	virtual void update_fields ( const Scalar dt ) override ;
	virtual void move_particles( const Scalar dt ) override ;

	virtual void adapt_meshes() override ;

	// Output
	virtual void dump_fields(unsigned frame) const override ;

private:
	PhaseMeshes  m_meshes ;
	std::unique_ptr<Phase>     m_grains ;

	DiphasicSolver m_solver ;

};


} //d6


#endif
