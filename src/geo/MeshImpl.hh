#ifndef D6_MESH_IMPL_HH
#define D6_MESH_IMPL_HH

#include "geo.fwd.hh"

#if( D6_MESH_IMPL == D6_MESH_TET_GRID )
#include "TetGrid.hh"
#elif( D6_MESH_IMPL == D6_MESH_OCTREE )
#include "Octree.hh"
#else
#include "Grid.hh"
#endif

#include "geo/MeshShapeFunction.hh"
#include "geo/P2ShapeFunction.hh"
#include "geo/UnstructuredShapeFunction.hh"

#endif
