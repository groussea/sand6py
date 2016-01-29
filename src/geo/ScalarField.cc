
#include "ScalarField.hh"
#include "FieldBase.impl.hh"

#include "Tensor.hh"

#include "Grid.hh"
#if HAS_TET
#include "TetGrid.hh"
#endif

#include "MeshShapeFunction.hh"

namespace d6
{

#if (D6_DIM==2)
template <typename ShapeFuncT>
void AbstractScalarField< ShapeFuncT >::get_spi_tensor(const Location &x, Mat &tensor) const
{
	VecR spi ; spi[0] = Base::eval_at(x) ;
	tensor_view( spi ).get( tensor ) ;
}

template <typename ShapeFuncT>
void AbstractScalarField< ShapeFuncT >::add_spi_tensor(const Location &x, Mat &tensor) const
{
	VecR spi ; spi[0] = Base::eval_at(x) ;
	tensor_view( spi ).add( tensor ) ;
}
#endif

template class FieldBase< AbstractScalarField< Linear<Grid> > > ;
template class AbstractScalarField< Linear<Grid> > ;
#if HAS_TET
template class FieldBase< AbstractScalarField< TetGrid > > ;
template class AbstractScalarField< TetGrid > ;
#endif
}
