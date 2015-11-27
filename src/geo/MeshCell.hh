#ifndef D6_MESH_CELL_HH
#define D6_MESH_CELL_HH

#include "geo.fwd.hh"

#if( D6_MESH_IMPL == D6_MESH_TET_GRID )
#include "Tet.hh"
#else
#include "Voxel.hh"
#endif

#endif
