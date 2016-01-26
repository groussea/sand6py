#include "MeshBase.hh"

#include "Grid.hh"
#if HAS_TET
#include "TetGrid.hh"
#endif

namespace d6 {

template <typename Derived>
Vec MeshBase< Derived>::pos( const Location& loc ) const
{
	CellGeo geo ;
	get_geo( loc.cell, geo ) ;
	return geo.pos( loc.coords ) ;
}

template class MeshBase< Grid > ;
#if HAS_TET
template class MeshBase< TetGrid > ;
#endif


}  //d6
