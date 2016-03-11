#ifndef D6_MONO_SIMU_HH
#define D6_MONO_SIMU_HH

#include "simu/Simu.hh"

#include "PhaseMeshes.hh"
#include "PhaseSolver.hh"

namespace d6 {

class MonoSimu : public Simu
{

public:

	MonoSimu( const Config& config, const char* base_dir ) ;
	virtual ~MonoSimu() ;

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
	// Useful for warm-starting stresses at frictional boundary conditions
	std::vector< RBStresses > m_rbStresses  ;

	PhaseSolver m_solver ;

};


} //d6


#endif
