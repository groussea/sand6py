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
template < typename ShapeFunc > class AbstractSkewTsField ;


class Grid ;
class TetGrid ;
#if( D6_MESH_IMPL == D6_MESH_TET_GRID )
typedef TetGrid MeshImpl ;
#else
typedef Grid    MeshImpl ;
#endif

template < typename MeshT > struct Linear ;
template < typename MeshT > struct DGLinear ;
template < typename MeshT > struct DGConstant ;
template < typename MeshT > struct P2 ;

typedef MeshBase< MeshImpl > MeshType ;

template < typename CellT, int Order >
struct QuadraturePoints ;

}

#endif
