
#include "ScalarField.hh"
#include "FieldBase.impl.hh"

#include "Tensor.hh"

#include "Grid.hh"
#if HAS_TET
#include "TetGrid.hh"
#endif

namespace d6
{

#if (D6_DIM==2)
template <typename MeshT>
void AbstractScalarField< MeshT >::get_spi_tensor(const Vec &x, Mat &tensor) const
{
	VecR spi ; spi[0] = Base::eval_at(x) ;
	tensor_view( spi ).get( tensor ) ;
}

template <typename MeshT>
void AbstractScalarField< MeshT >::add_spi_tensor(const Vec &x, Mat &tensor) const
{
	VecR spi ; spi[0] = Base::eval_at(x) ;
	tensor_view( spi ).add( tensor ) ;
}
#endif

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
#if HAS_TET
template class FieldBase< AbstractScalarField< TetGrid > > ;
template class AbstractScalarField< TetGrid > ;
#endif
}
