#ifndef D6_MESH_TRAITS_HH
#define D6_MESH_TRAITS_HH


#include "geo/Voxel.hh"
#include "geo/Tet.hh"

#include <vector>

namespace d6 {

template <typename Derived >
struct MeshTraits {
};

// Grid

class  Grid ;
struct GridIterator ;

template < >
struct MeshTraits< Grid > {
	typedef GridIterator CellIterator ;

	typedef ArrWi Cell    ;
	typedef Voxel CellGeo ;
	static constexpr Index NV = CellGeo::NV ;
	static constexpr Index NC = CellGeo::NC ;
	static constexpr Index NQ = CellGeo::NQ ;

	typedef std::vector<Cell> Cells ;
};

// TetGrid

class  TetGrid ;
struct TetGridIterator ;

template < >
struct MeshTraits< TetGrid > {
	typedef TetGridIterator CellIterator ;

	typedef Eigen::Array< int, WD+1, 1 > Cell  ;
	typedef Tet CellGeo ;
	static constexpr Index NV = CellGeo::NV ;
	static constexpr Index NC = CellGeo::NC ;
	static constexpr Index NQ = CellGeo::NQ ;

	typedef std::vector<Cell> Cells ;
};

} //d6

#endif
