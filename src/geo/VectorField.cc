
#include "VectorField.hh"
#include "FieldBase.impl.hh"

#include "Tensor.hh"

#include "Grid.hh"
#include "TetGrid.hh"

#include "MeshShapeFunction.hh"

namespace d6
{

template class FieldBase< AbstractVectorField< Linear<Grid> > > ;
template class AbstractVectorField< Linear<Grid> > ;
template class FieldBase< AbstractVectorField< Linear< TetGrid > > > ;
template class AbstractVectorField< Linear< TetGrid > > ;

template class FieldBase< AbstractVectorField< DGLinear<Grid> > > ;
template class AbstractVectorField< DGLinear<Grid> > ;
template class FieldBase< AbstractVectorField< DGLinear< TetGrid > > > ;
template class AbstractVectorField< DGLinear< TetGrid > > ;

}

