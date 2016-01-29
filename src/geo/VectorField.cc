
#include "VectorField.hh"
#include "FieldBase.impl.hh"

#include "Tensor.hh"

#include "Grid.hh"
#if HAS_TET
#include "TetGrid.hh"
#endif

#include "MeshShapeFunction.hh"

namespace d6
{

#if (D6_DIM==3)
template <typename ShapeFuncT>
void AbstractVectorField< ShapeFuncT >::get_spi_tensor(const Vec &x, Mat &tensor) const
{
	tensor_view( Base::eval_at(x) ).get( tensor ) ;
}

template <typename ShapeFuncT>
void AbstractVectorField< ShapeFuncT >::add_spi_tensor(const Vec &x, Mat &tensor) const
{
	tensor_view( Base::eval_at(x) ).add( tensor ) ;
}
#endif

template class FieldBase< AbstractVectorField< Linear<Grid> > > ;
template class AbstractVectorField< Linear<Grid> > ;
#if HAS_TET
template class FieldBase< AbstractVectorField< Linear< TetGrid > > > ;
template class AbstractVectorField< Linear< TetGrid > > ;
#endif

}

