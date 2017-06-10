#include "FluidPhase.hh"


#include "geo/FieldBase.impl.hh" //interpolate

namespace d6 {


FluidPhase::FluidPhase( const PhaseMeshes & meshes, const FluidPhase& src )
	: pressure(src.pressure.interpolate<PrimalShape>(meshes.primal())),
	  velocity(src.velocity.interpolate<PrimalShape>(meshes.primal())),
	  mavg_vel(src.mavg_vel.interpolate<PrimalShape>(meshes.primal())),
	  sym_grad(src.sym_grad.interpolate<PrimalShape>(meshes.primal()))

{
}


} //d6
