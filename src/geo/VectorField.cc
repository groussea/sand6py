
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

template class FieldBase< AbstractVectorField< Linear<Grid> > > ;
template class AbstractVectorField< Linear<Grid> > ;
#if HAS_TET
template class FieldBase< AbstractVectorField< Linear< TetGrid > > > ;
template class AbstractVectorField< Linear< TetGrid > > ;
#endif

}

