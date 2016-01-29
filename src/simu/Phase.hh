#ifndef D6_PHASE_DESCR_HH
#define D6_PHASE_DESCR_HH

#include "geo/geo.fwd.hh"
#include "geo/ScalarField.hh"
#include "geo/VectorField.hh"
#include "geo/TensorField.hh"

#include "geo/MeshImpl.hh"

namespace d6 {

struct Phase
{

	ScalarField fraction ;
	VectorField velocity ;

	TensorField stresses ;
	TensorField sym_grad ;
#if D6_DIM == 3
	VectorField spi_grad ;
#else
	ScalarField spi_grad ;
#endif

	VectorField grad_phi ;
	VectorField geo_proj ;
	VectorField fcontact ;

	Phase( const MeshType& mesh )
		: fraction(Linear<MeshImpl>(mesh)), velocity(Linear<MeshImpl>(mesh)),
		  stresses(Linear<MeshImpl>(mesh)), sym_grad(Linear<MeshImpl>(mesh)),
		  spi_grad(Linear<MeshImpl>(mesh)), grad_phi(Linear<MeshImpl>(mesh)),
		  geo_proj(Linear<MeshImpl>(mesh)), fcontact(Linear<MeshImpl>(mesh))
	{}

	template < typename Archive >
	void serialize( Archive &ar, unsigned int ) {
		ar & fraction ;
		ar & velocity ;
		ar & stresses ;
		ar & sym_grad ;
		ar & spi_grad ;
		ar & grad_phi ;
		ar & fcontact ;
	}
};

} //d6

#endif
