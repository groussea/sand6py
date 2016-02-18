#ifndef D6_PHASE_DESCR_HH
#define D6_PHASE_DESCR_HH

#include "simu/PhaseFields.hh"

#include "geo/ScalarField.hh"
#include "geo/VectorField.hh"
#include "geo/TensorField.hh"

#include "geo/MeshImpl.hh"

namespace d6 {

struct Phase
{


	PrimalScalarField fraction ;
	PrimalVectorField velocity ;

	DualTensorField stresses ;
	DualTensorField sym_grad ;
	DualSkewTsField spi_grad ;

	PrimalVectorField grad_phi ;
	PrimalVectorField geo_proj ;
	PrimalVectorField fcontact ;

	Phase( const PhaseMeshes & meshes )
		: fraction(meshes.primal()), velocity(meshes.primal()),
		  stresses(meshes.  dual()), sym_grad(meshes.  dual()),
		  spi_grad(meshes.  dual()), grad_phi(meshes.primal()),
		  geo_proj(meshes.primal()), fcontact(meshes.primal())
	{}

	Phase( const PhaseMeshes & meshes, const Phase& src ) ;

	template < typename Archive >
	void serialize( Archive &ar, unsigned int ) {
		ar & fraction ;
		ar & velocity ;
		ar & stresses ;
		ar & sym_grad ;
		ar & spi_grad ;
		ar & grad_phi ;
		ar & fcontact ;
		ar & geo_proj ;
	}
};

} //d6

#endif
