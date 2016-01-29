#ifndef D6_GEO_FWD_HH
#define D6_GEO_FWD_HH

#define D6_MESH_GRID     0
#define D6_MESH_TET_GRID 1

#ifndef D6_MESH_IMPL
#define D6_MESH_IMPL D6_MESH_GRID
#endif

namespace d6 {

template < typename M > class MeshBase ;
template < typename S > class ShapeFuncBase ;


template< typename ValueType > struct Expr ;

template< typename Derived > class FieldBase ;
template< typename Derived > struct FieldTraits ;
template < typename ShapeFunc > class AbstractScalarField ;
template < typename ShapeFunc > class AbstractVectorField ;
template < typename ShapeFunc > class AbstractTensorField ;


class Grid ;
class TetGrid ;
#if( D6_MESH_IMPL == D6_MESH_TET_GRID )
typedef TetGrid MeshImpl ;
#else
typedef Grid    MeshImpl ;
#endif

template < typename MeshT > struct Linear ;
template < typename MeshT > struct DGLinear ;

typedef AbstractScalarField< Linear<MeshImpl> >  ScalarField ;
typedef AbstractVectorField< Linear<MeshImpl> >  VectorField ;
typedef AbstractTensorField< Linear<MeshImpl> >  TensorField ;


typedef MeshBase< MeshImpl > MeshType ;

}

#endif
