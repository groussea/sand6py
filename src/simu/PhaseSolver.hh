#ifndef D6_PHASE_SOLVE_HH
#define D6_PHASE_SOLVE_HH

#include "geo/BoundaryInfo.hh"
#include "geo/MeshBase.hh"

#include "utils/scalar.hh"

#include <vector>

namespace d6 {

class DynParticles ;
class Config ;
struct Phase ;

class PhaseSolver {

public:
	explicit PhaseSolver( const DynParticles& particles ) ;

	void step(const Config &config, Phase& phase ) ;

private:

	void computeActiveNodes(const MeshType &mesh,
							const std::vector< bool >& activeCells ) ;

	struct Active {

		static const Index s_Inactive  ;

		Index nNodes ;
		typename MeshType::Cells cells ;
		std::vector< Index > indices ;

		Active()	 : nNodes( 0 ) {}

		void reset( Index totNodes )
		{
			nNodes = 0 ;
			cells.clear();
			indices.assign( totNodes, s_Inactive );
		}

		Index count() const { return nNodes ; }
	};


	const DynParticles& m_particles ;

	Active m_phaseNodes ;
	Active m_couplingNodes ;

	BoundaryConditions m_surfaceNodes ;
};


} //d6


#endif
