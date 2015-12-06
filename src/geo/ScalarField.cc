
#include "ScalarField.hh"
#include "FieldBase.impl.hh"

#include "Grid.hh"
#include "TetGrid.hh"

namespace d6
{

template <typename MeshT >
Vec AbstractScalarField< MeshT >::grad_at( const Vec& x ) const 
{
	typename MeshType::Location loc ;
	typename MeshType::NodeList nodes ;
	typename MeshType::Derivatives dc_dx ;

	m_mesh.locate( x, loc );
	m_mesh.list_nodes( loc.cell, nodes );
	m_mesh.get_derivatives( loc, dc_dx ) ;

	// v(x) = sum c_k(x) v_k

	Vec grad = Vec::Zero() ;
	for( Index k = 0 ; k < nodes.rows() ; ++k ) {
		grad += dc_dx.row(k) * Base::segment( nodes[k] ) ;
	}
	
	return grad ;
}

template class FieldBase< AbstractScalarField< Grid > > ;
template class AbstractScalarField< Grid > ;
template class FieldBase< AbstractScalarField< TetGrid > > ;
template class AbstractScalarField< TetGrid > ;

}
