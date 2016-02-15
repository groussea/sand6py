
#include "VectorField.hh"
#include "FieldBase.impl.hh"

#include "instanciations.hh"

namespace d6
{

#define INSTANTIATE( Shape ) \
	template class FieldBase< AbstractVectorField< Shape > > ; \
	template class AbstractVectorField< Shape > ;

EXPAND_INSTANTIATIONS
#undef INSTANTIATE

}

