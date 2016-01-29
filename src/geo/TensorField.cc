
#include "TensorField.hh"
#include "FieldBase.impl.hh"

#include "Tensor.hh"

#include "Grid.hh"
#if HAS_TET
#include "TetGrid.hh"
#endif

#include "MeshShapeFunction.hh"

namespace d6
{

template <typename ShapeFuncT>
void AbstractTensorField< ShapeFuncT >::get_sym_tensor(const Location &x, Mat &tensor) const
{
	tensor_view( Base::eval_at(x) ).get( tensor ) ;
}

template <typename ShapeFuncT>
void AbstractTensorField< ShapeFuncT >::add_sym_tensor(const Location &x, Mat &tensor) const
{
	tensor_view( Base::eval_at(x) ).add( tensor ) ;
}

template class FieldBase< AbstractTensorField< Linear<Grid> > > ;
template class AbstractTensorField< Linear<Grid> > ;
#if HAS_TET
template class FieldBase< AbstractTensorField< Linear<TetGrid> > > ;
template class AbstractTensorField< Linear<TetGrid> > ;
#endif
}

