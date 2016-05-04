#include "FluidPhase.hh"


#include "geo/FieldBase.impl.hh" //interpolate

namespace d6 {


FluidPhase::FluidPhase( const PhaseMeshes & meshes, const FluidPhase& src )
	: pressure(src.pressure.interpolate<PrimalShape>(meshes.primal())),
	  velocity(src.velocity.interpolate<PrimalShape>(meshes.primal()))

{
}


} //d6
