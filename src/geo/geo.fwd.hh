#ifndef D6_GEO_FWD_HH
#define D6_GEO_FWD_HH

#define D6_MESH_GRID     0
#define D6_MESH_TET_GRID 1

#ifndef D6_MESH_IMPL
#define D6_MESH_IMPL D6_MESH_GRID
#endif

namespace d6 {

class Grid ;
template < typename M > class MeshBase ;

template< typename ValueType > struct Expr ;

template< typename Derived > class FieldBase ;
template< typename Derived > struct FieldTraits ;
template < typename MeshT > class AbstractScalarField ;
template < typename MeshT > class AbstractVectorField ;
template < typename MeshT > class AbstractTensorField ;


#if( D6_MESH_IMPL == D6_MESH_TET_GRID )
typedef TetGrid MeshImpl ;
#else
typedef Grid    MeshImpl ;
#endif



typedef MeshBase< MeshImpl > MeshType ;

typedef AbstractScalarField< MeshImpl >  ScalarField ;
typedef AbstractVectorField< MeshImpl >  VectorField ;
typedef AbstractTensorField< MeshImpl >  TensorField ;


}

#endif
