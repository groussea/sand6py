#ifndef D6_PHASE_FIELDS_HH
#define D6_PHASE_FIELDS_HH

#include "geo/geo.fwd.hh"

namespace d6 {

    typedef   Linear<MeshImpl> PrimalShape ;
    typedef DGLinear<MeshImpl>   DualShape ;

    typedef AbstractScalarField< PrimalShape > PrimalScalarField ;
    typedef AbstractVectorField< PrimalShape > PrimalVectorField ;
    typedef AbstractTensorField< PrimalShape > PrimalTensorField ;

    typedef AbstractScalarField<   DualShape > DualScalarField ;
    typedef AbstractTensorField<   DualShape > DualTensorField ;
    typedef AbstractSkewTsField<   DualShape > DualSkewTsField ;

    typedef PrimalTensorField RBStresses ;


} //d6


#endif
