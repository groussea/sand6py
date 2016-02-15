
#include "ScalarField.hh"
#include "FieldBase.impl.hh"

#include "geo/instanciations.hh"

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

#define INSTANTIATE( Shape ) \
	template class FieldBase< AbstractScalarField< Shape > > ; \
	template class AbstractScalarField< Shape > ;

EXPAND_INSTANTIATIONS
#undef INSTANTIATE

}
