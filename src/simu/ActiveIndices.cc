#include "ActiveIndices.hh"

#include "PhaseFields.hh"

#include "geo/ScalarField.hh"
#include "geo/VectorField.hh"
#include "geo/TensorField.hh"

#include "geo/MeshImpl.hh"

namespace d6 {

const Index Active::s_Inactive = -1 ;

void Active::computeRevIndices()
{
	revIndices.resize( nNodes );

#pragma omp parallel for
	for( size_t i = 0 ; i < indices.size() ; ++i  ) {
		const Index idx = indices[ i ] ;
		if( idx != Active::s_Inactive ) {
			revIndices[idx-offset] = i ;
		}
	}
}

void Active::setOffset(const Index o)
{
	offset = o ;

#pragma omp parallel for
	for( size_t i = 0 ; i < indices.size() ; ++i  ) {
		if( indices[i] != Active::s_Inactive )
			indices[i] += o ;
	}
}

} //d6
