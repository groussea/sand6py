#ifndef D6_MESH_TRAITS_HH
#define D6_MESH_TRAITS_HH

#define HAS_TET (D6_DIM == 3)

#include "geo/Voxel.hh"
#if HAS_TET
#include "geo/Tet.hh"
#endif

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

#if HAS_TET
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
#endif

} //d6

#endif
