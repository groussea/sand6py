
#include "TensorField.hh"
#include "FieldBase.impl.hh"

#include "Tensor.hh"

#include "Grid.hh"
#if HAS_TET
#include "TetGrid.hh"
#endif

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
#if HAS_TET
template class FieldBase< AbstractTensorField< TetGrid > > ;
template class AbstractTensorField< TetGrid > ;
#endif
}

