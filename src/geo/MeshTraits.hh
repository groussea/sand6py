#ifndef D6_MESH_TRAITS_HH
#define D6_MESH_TRAITS_HH

#include "Voxel.hh"
#include <vector>

namespace d6 {

template <typename Derived >
struct MeshTraits {
};

class  Grid ;
struct GridIterator ;
template < >
struct MeshTraits< Grid > {
	typedef GridIterator CellIterator ;

	typedef Vec3i Cell    ;
	typedef Voxel CellGeo ;
	static constexpr Index NV = CellGeo::NV ;
	static constexpr Index NC = CellGeo::NC ;
	static constexpr Index NQ = CellGeo::NQ ;

	typedef std::vector<Cell> Cells ;
};

} //d6

#endif
