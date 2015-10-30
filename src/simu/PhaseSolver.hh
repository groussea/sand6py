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
struct PhaseMatrices ;

class PhaseSolver {

public:
	explicit PhaseSolver( const DynParticles& particles ) ;

	void step(const Config &config, Phase& phase ) ;

private:

	void computeActiveNodes(const MeshType &mesh,
							const std::vector< bool >& activeCells ) ;

	void computeProjectors( PhaseMatrices& matrices ) const ;

	void assembleMatrices( const Config& c, const MeshType& mesh,
						   const DynVec &phiInt,
						   PhaseMatrices& matrices ) const ;

	void solveComplementarity(const Config&c, const PhaseMatrices& matrices , const DynVec &fraction,
							  DynVec &u, Phase &phase) const ;

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
		Index origSize() const { return indices.size() ; }

		template < typename Derived >
		void field2var( const FieldBase<Derived> &field, DynVec & var ) const ;
		template < typename Derived >
		void var2field( const DynVec & var, FieldBase<Derived> &field ) const ;
	};

	const DynParticles& m_particles ;

	Active m_phaseNodes ;
	Active m_couplingNodes ;

	BoundaryConditions m_surfaceNodes ;
};


} //d6


#endif
