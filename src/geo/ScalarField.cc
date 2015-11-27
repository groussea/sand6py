
#include "ScalarField.hh"
#include "FieldBase.impl.hh"

#include "Grid.hh"
#include "TetGrid.hh"

namespace d6
{

template class FieldBase< AbstractScalarField< Grid > > ;
template class AbstractScalarField< Grid > ;
template class FieldBase< AbstractScalarField< TetGrid > > ;
template class AbstractScalarField< TetGrid > ;

}
