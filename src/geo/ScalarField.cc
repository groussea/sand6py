
#include "ScalarField.hh"
#include "FieldBase.impl.hh"

#include "Tensor.hh"

#include "Grid.hh"
#include "TetGrid.hh"

#include "MeshShapeFunction.hh"
#include "P2ShapeFunction.hh"

namespace d6
{


template<typename ShapeFuncT >
Vec AbstractScalarField< ShapeFuncT >::grad_at( const typename ShapeFuncType::Location& loc ) const
{
	typename ShapeFuncType::NodeList nodes ;
	typename ShapeFuncType::Derivatives dc_dx ;

	Base::shape().list_nodes( loc, nodes );
	Base::shape().get_derivatives( loc, dc_dx ) ;

	// v(x) = sum c_k(x) v_k

	Vec grad = Vec::Zero() ;
	for( Index k = 0 ; k < nodes.rows() ; ++k ) {
		grad += dc_dx.row(k) * Base::segment( nodes[k] ) ;
	}

	return grad ;
}


template class FieldBase< AbstractScalarField< Linear<Grid> > > ;
template class AbstractScalarField< Linear<Grid> > ;
template class FieldBase< AbstractScalarField< Linear<TetGrid> > > ;
template class AbstractScalarField< Linear<TetGrid> > ;


template class FieldBase< AbstractScalarField< DGLinear<Grid> > > ;
template class AbstractScalarField< DGLinear<Grid> > ;
template class FieldBase< AbstractScalarField< DGLinear<TetGrid> > > ;
template class AbstractScalarField< DGLinear<TetGrid> > ;
template class FieldBase< AbstractScalarField< P2<TetGrid> > > ;
template class AbstractScalarField< P2<TetGrid> > ;

}
