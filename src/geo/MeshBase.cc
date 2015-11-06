#include "MeshBase.hh"

#include "Grid.hh"

namespace d6 {

template <typename Derived>
Vec MeshBase< Derived>::pos( const Location& loc ) const
{
	CellGeo geo ;
	get_geo( loc.cell, geo ) ;
	return geo.pos( loc.coords ) ;
}

template class MeshBase< Grid > ;


}  //d6
