
#include "VectorField.hh"
#include "FieldBase.impl.hh"

#include "instanciations.hh"

namespace d6
{

template<typename ShapeFuncT >
Mat AbstractVectorField< ShapeFuncT >::grad_at( const typename ShapeFuncType::Location& loc ) const
{
	typename ShapeFuncType::NodeList nodes ;
	typename ShapeFuncType::Derivatives dc_dx ;

	Base::shape().list_nodes( loc, nodes );
	Base::shape().get_derivatives( loc, dc_dx ) ;

	// v(x) = sum c_k(x) v_k

	Mat grad = Mat::Zero() ;
	for( Index k = 0 ; k < nodes.rows() ; ++k ) {
		grad += Base::segment( nodes[k] ) * dc_dx.row(k) ;
	}

	return grad ;
}


#define INSTANTIATE( Shape ) \
	template class FieldBase< AbstractVectorField< Shape > > ; \
	template class AbstractVectorField< Shape > ;

EXPAND_INSTANTIATIONS
#undef INSTANTIATE

}

