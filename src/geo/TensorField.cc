
#include "TensorField.hh"
#include "FieldBase.impl.hh"

#include "Grid.hh"
#include "Tensor.hh"

namespace d6
{

template <typename MeshT>
void AbstractTensorField< MeshT >::get_sym_tensor(const Vec &x, Mat &tensor) const
{
	tensor_view( Base::eval_at(x) ).get( tensor ) ;
}

template <typename MeshT>
void AbstractTensorField< MeshT >::add_sym_tensor(const Vec &x, Mat &tensor) const
{
	tensor_view( Base::eval_at(x) ).add( tensor ) ;
}

template class FieldBase< AbstractTensorField< Grid > > ;
template class AbstractTensorField< Grid > ;

}

