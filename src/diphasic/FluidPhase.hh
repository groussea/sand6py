#ifndef D6_FLUID_PHASE_HH
#define D6_FLUID_PHASE_HH


#include "mono/PhaseMeshes.hh"

#include "geo/ScalarField.hh"
#include "geo/VectorField.hh"
#include "geo/TensorField.hh"

#include "geo/MeshImpl.hh"

namespace d6 {

struct FluidPhase
{
	PrimalScalarField pressure ;
	PrimalVectorField velocity ;

	FluidPhase( const PhaseMeshes & meshes )
		: pressure(meshes.primal()), velocity(meshes.primal())
	{}

	FluidPhase( const PhaseMeshes & meshes, const FluidPhase& src ) ;

	template < typename Archive >
	void serialize( Archive &ar, unsigned int ) {
		ar & pressure ;
		ar & velocity ;
	}

};

} //d6

#endif
