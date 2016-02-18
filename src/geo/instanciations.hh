#ifndef D6_INSTANTIATIONS_HH
#define D6_INSTANTIATIONS_HH

#include "Grid.hh"
#include "Octree.hh"
#include "TetGrid.hh"

#include "MeshShapeFunction.hh"
#include "P2ShapeFunction.hh"
#include "UnstructuredShapeFunction.hh"


#define EXPAND_INSTANTIATIONS \
	INSTANTIATE(     Linear<Grid   >  ) \
	INSTANTIATE(   DGLinear<Grid   >  ) \
	INSTANTIATE( DGConstant<Grid   >  ) \
	INSTANTIATE(     Linear<Octree >  ) \
	INSTANTIATE(   DGLinear<Octree >  ) \
	INSTANTIATE( DGConstant<Octree >  ) \
	INSTANTIATE(     Linear<TetGrid>  ) \
	INSTANTIATE(   DGLinear<TetGrid>  ) \
	INSTANTIATE( DGConstant<TetGrid>  ) \
	INSTANTIATE( UnstructuredShapeFunc  ) \
	EXPAND_INSTANTIATIONS_DIM

#define EXPAND_INSTANTIATIONS_2D \
	INSTANTIATE( P2<TetGrid> ) \

#define EXPAND_INSTANTIATIONS_3D \


#if D6_DIM == 2
	#define EXPAND_INSTANTIATIONS_DIM EXPAND_INSTANTIATIONS_2D
#else
	#define EXPAND_INSTANTIATIONS_DIM EXPAND_INSTANTIATIONS_3D
#endif

#endif
