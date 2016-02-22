#ifndef D6_PHASE_FIELDS_HH
#define D6_PHASE_FIELDS_HH

#include "geo/geo.fwd.hh"

#include <memory>

namespace d6 {

typedef MeshImpl PrimalMesh ;

typedef Linear<PrimalMesh> PrimalShape ;
//typedef   P2<PrimalMesh> PrimalShape ;

#ifdef D6_UNSTRUCTURED_DUAL
typedef UnstructuredDOFs      DualMesh ;
typedef UnstructuredShapeFunc DualShape ;
#else
typedef MeshImpl   DualMesh ;

#ifdef D6_DG_STRESSES
typedef DGLinear<  DualMesh> DualShape ;
//typedef DGConstant<  DualMesh> DualShape ;
#else
typedef   Linear<  DualMesh> DualShape ;
#endif
#endif

typedef AbstractScalarField< PrimalShape > PrimalScalarField ;
typedef AbstractVectorField< PrimalShape > PrimalVectorField ;
typedef AbstractTensorField< PrimalShape > PrimalTensorField ;

typedef AbstractScalarField<   DualShape > DualScalarField ;
typedef AbstractTensorField<   DualShape > DualTensorField ;
typedef AbstractSkewTsField<   DualShape > DualSkewTsField ;

typedef PrimalTensorField RBStresses ;



} //d6


#endif
