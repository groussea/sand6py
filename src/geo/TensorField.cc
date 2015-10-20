
#include "TensorField.hh"
#include "FieldBase.impl.hh"

#include "Grid.hh"

namespace d6
{

template class FieldBase< AbstractTensorField< Grid > > ;
template class AbstractTensorField< Grid > ;

}

