
#include "utils/Config.hh"
#include "geo/ScalarField.hh"

#include "geo/MeshShapeFunction.hh"
#include "geo/UnstructuredShapeFunction.hh"
#include "geo/Grid.hh"

#include "geo/FieldBase.impl.hh"

#include <gtest/gtest.h>

using namespace d6 ;

TEST( geo, interp )
{
	Config c ;

	Grid g( c.box, c.res ) ;

	AbstractScalarField< Linear<Grid> > lin ( g ) ;
	lin.set_constant( 1 ) ;


	AbstractScalarField< DGLinear<Grid> > dg ( lin.interpolate< DGLinear<Grid> >() ) ;

	ASSERT_EQ( g.nCells() * Voxel::NQ, dg.size() ) ;
	ASSERT_DOUBLE_EQ( 1., dg( c.box/2 ) ) ;
	ASSERT_DOUBLE_EQ( dg.size(), dg.flatten().sum() ) ;

	lin = dg.interpolate< Linear<Grid> > () ;

	ASSERT_DOUBLE_EQ( lin.size(), lin.flatten().sum() ) ;

	lin = lin.interpolate< Linear<Grid> > () ;

}

TEST(geo, shape_functions)
{
	typename UnstructuredDOFs::Vertices vertices ;
	typename UnstructuredDOFs::Weights weights ;

	UnstructuredDOFs dofs ( vertices, weights) ;
	UnstructuredShapeFunc shape( dofs ) ;
	AbstractScalarField< UnstructuredShapeFunc > prtScalarField( shape ) ;


}
