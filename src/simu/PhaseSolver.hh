#ifndef D6_PHASE_SOLVE_HH
#define D6_PHASE_SOLVE_HH

#include "ActiveIndices.hh"

#include "geo/BoundaryInfo.hh"
#include "geo/MeshBase.hh"

#include "utils/scalar.hh"

#include <vector>

namespace d6 {

class DynParticles ;
struct Phase ;
struct PhaseMatrices ;

class RigidBody ;
struct RigidBodyData ;

struct Config ;
class Stats ;

class PhaseSolver {

public:
	explicit PhaseSolver(
			const DynParticles& particles
			) ;

	~PhaseSolver() ;

	void step(const Config &config, Phase& phase, Stats &stats,
			  std::vector<RigidBody> &rigidBodies,
			  std::vector<TensorField> &rbStresses
			  ) ;

private:

	void computeActiveNodes(const std::vector< bool >& activeCells,
							const VectorField &grad_phi ) ;
	void computeActiveBodies( std::vector<RigidBody> &rigidBodies,
							  std::vector<TensorField> &rbStresses,
							  std::vector< RigidBodyData >& rbData ) ;



	void computeProjectors( PhaseMatrices& matrices, const bool weakStressBC,
							const std::vector< RigidBodyData >& rbData
							) const ;

	void computeAnisotropy(const DynVec& orientation,
						   const Config &config, PhaseMatrices& matrices) const ;

	void assembleMatrices( const Config& c, const MeshType& mesh,
						   const ScalarField &phiInt,
						   PhaseMatrices& matrices,
						   std::vector< RigidBodyData >& rbData
						   ) const ;

	void solveComplementarity(const Config&c, const PhaseMatrices& matrices ,
							  std::vector< RigidBodyData >& rbData,
							  const DynVec &fraction,
							  const DynVec &cohesion, const DynVec &inertia,
							  DynVec &u, Phase &phase, Stats &stats ) const ;
	void enforceMaxFrac(const Config &c, const PhaseMatrices &matrices,
									   std::vector<RigidBodyData> &rbData,
									   const DynVec &fraction,
									   DynVec &depl ) const ;

	void computeGradPhi( const ScalarField& fraction, const ScalarField& volumes, VectorField &grad_phi ) const ;


	Index nSuppNodes() const
	{
		return m_totRbNodes ;
	}


	const DynParticles& m_particles ;

	Active m_phaseNodes ;
	Active m_couplingNodes ;

	BoundaryConditions m_surfaceNodes ;

	Index m_totRbNodes ;
};


} //d6


#endif
