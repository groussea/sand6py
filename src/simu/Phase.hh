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
	VectorField spi_grad ;

	VectorField grad_phi ;
	VectorField geo_proj ;
	VectorField fcontact ;

	Phase( const MeshType& mesh )
		: fraction(mesh), velocity(mesh),
		  stresses(mesh), sym_grad(mesh),
		  spi_grad(mesh), grad_phi(mesh),
		  geo_proj(mesh), fcontact(mesh)
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
