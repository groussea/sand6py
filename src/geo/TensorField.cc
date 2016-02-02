
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

template class FieldBase< AbstractTensorField< Linear<Grid> > > ;
template class AbstractTensorField< Linear<Grid> > ;
template class FieldBase< AbstractSkewTsField< Linear<Grid> > > ;
template class AbstractSkewTsField< Linear<Grid> > ;

#if HAS_TET
template class FieldBase< AbstractTensorField< Linear<TetGrid> > > ;
template class AbstractTensorField< Linear<TetGrid> > ;
template class FieldBase< AbstractSkewTsField< Linear<TetGrid> > > ;
template class AbstractSkewTsField< Linear<TetGrid> > ;
#endif
}

