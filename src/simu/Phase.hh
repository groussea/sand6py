#ifndef D6_PHASE_DESCR_HH
#define D6_PHASE_DESCR_HH

#include "geo/geo.fwd.hh"
#include "geo/ScalarField.hh"
#include "geo/VectorField.hh"
#include "geo/TensorField.hh"

namespace d6 {

struct Phase
{

	ScalarField fraction ;
	VectorField velocity ;

	TensorField stresses ;
	TensorField sym_grad ;
	VectorField spi_grad ;

	Phase( const MeshType& mesh )
		: fraction(mesh), velocity(mesh),
		  stresses(mesh), sym_grad(mesh),
		  spi_grad(mesh)
	{}

};

} //d6

#endif
