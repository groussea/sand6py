
#include "VectorField.hh"
#include "FieldBase.impl.hh"

#include "Tensor.hh"

#include "Grid.hh"
#include "TetGrid.hh"

namespace d6
{

template <typename MeshT>
void AbstractVectorField< MeshT >::get_spi_tensor(const Vec &x, Mat &tensor) const
{
	tensor_view( Base::eval_at(x) ).get( tensor ) ;
}

template <typename MeshT>
void AbstractVectorField< MeshT >::add_spi_tensor(const Vec &x, Mat &tensor) const
{
	tensor_view( Base::eval_at(x) ).add( tensor ) ;
}

template class FieldBase< AbstractVectorField< Grid > > ;
template class AbstractVectorField< Grid > ;
template class FieldBase< AbstractVectorField< TetGrid > > ;
template class AbstractVectorField< TetGrid > ;

}

