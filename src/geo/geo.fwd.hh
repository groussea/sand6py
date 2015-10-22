#ifndef D6_GEO_FWD_HH
#define D6_GEO_FWD_HH

namespace d6 {

class Grid ;
template < typename M > class MeshBase ;

template< typename Derived > struct FieldBase ;
template < typename MeshT > class AbstractScalarField ;
template < typename MeshT > class AbstractVectorField ;
template < typename MeshT > class AbstractTensorField ;

typedef Grid MeshImpl ;
typedef MeshBase< MeshImpl > MeshType ;

typedef AbstractScalarField< MeshImpl >  ScalarField ;
typedef AbstractVectorField< MeshImpl >  VectorField ;
typedef AbstractTensorField< MeshImpl >  TensorField ;


}

#endif
