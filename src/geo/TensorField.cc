
#include "TensorField.hh"
#include "FieldBase.impl.hh"

#include "Tensor.hh"

#include "instanciations.hh"

namespace d6
{

template <typename ShapeFuncT>
void AbstractTensorField< ShapeFuncT >::get_sym_tensor(const Location &x, Mat &tensor) const
{
	const typename Base::ValueType& v = Base::eval_at(x) ;
	tensor_view( v ).get( tensor ) ;
}

template <typename ShapeFuncT>
void AbstractTensorField< ShapeFuncT >::add_sym_tensor(const Location &x, Mat &tensor) const
{
	const typename Base::ValueType& v = Base::eval_at(x) ;
	tensor_view( v ).add( tensor ) ;
}

template <typename ShapeFuncT>
void AbstractSkewTsField< ShapeFuncT >::get_spi_tensor(const Location &x, Mat &tensor) const
{
	const typename Base::ValueType& v = Base::eval_at(x) ;
	tensor_view( v ).get( tensor ) ;
}

template <typename ShapeFuncT>
void AbstractSkewTsField< ShapeFuncT >::add_spi_tensor(const Location &x, Mat &tensor) const
{
	const typename Base::ValueType& v = Base::eval_at(x) ;
	tensor_view( v ).add( tensor ) ;
}

#define INSTANTIATE( Shape ) \
	template class FieldBase< AbstractTensorField< Shape > > ; \
	template class AbstractTensorField< Shape > ; \
	template class FieldBase< AbstractSkewTsField< Shape > > ; \
	template class AbstractSkewTsField< Shape > ;

EXPAND_INSTANTIATIONS
#undef INSTANTIATE

}

