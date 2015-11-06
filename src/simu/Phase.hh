#ifndef D6_PHASE_DESCR_HH
#define D6_PHASE_DESCR_HH

#include "geo/geo.fwd.hh"
#include "geo/ScalarField.hh"
#include "geo/VectorField.hh"
#include "geo/TensorField.hh"

#include "geo/Grid.hh"

namespace d6 {

struct Phase
{

	ScalarField fraction ;
	VectorField velocity ;

	TensorField stresses ;
	TensorField sym_grad ;
	VectorField spi_grad ;

	VectorField geo_proj ;

	Phase( const MeshType& mesh )
		: fraction(mesh), velocity(mesh),
		  stresses(mesh), sym_grad(mesh),
		  spi_grad(mesh), geo_proj(mesh)
	{}

	template < typename Archive >
	void serialize( Archive &ar, unsigned int ) {
		ar & fraction ;
		ar & velocity ;
		ar & stresses ;
		ar & sym_grad ;
		ar & spi_grad ;
	}
};

} //d6

#endif
