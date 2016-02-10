
#include "utils/Config.hh"
#include "geo/ScalarField.hh"

#include "geo/MeshShapeFunction.hh"
#include "geo/UnstructuredShapeFunction.hh"
#include "geo/Grid.hh"

#include "geo/FieldBase.impl.hh"

#include <gtest/gtest.h>

using namespace d6 ;

TEST( geo, rvalue )
{
	Config c ;
	Grid g( c.box, c.res ) ;

	AbstractScalarField< Linear<Grid> > field(g) ;

	const Scalar* orig_ptr = field.flatten().data() ;

	AbstractScalarField< Linear<Grid> > field_copy( field ) ;
	const Scalar* copy_ptr = field_copy.flatten().data() ;

	ASSERT_NE( orig_ptr, copy_ptr ) ;
	field_copy.flatten().swap( field.flatten() ) ;

	copy_ptr = field_copy.flatten().data() ;
	ASSERT_EQ( orig_ptr, copy_ptr ) ;

	AbstractScalarField< Linear<Grid> > field_move( std::move(field_copy) ) ;
	const Scalar* move_ptr = field_move.flatten().data() ;

	ASSERT_EQ( orig_ptr, move_ptr ) ;


}

TEST( geo, interp )
{
	Config c ;

	Grid g( c.box, c.res ) ;

	AbstractScalarField< Linear<Grid> > lin ( g ) ;
	lin.set_constant( 1 ) ;


	AbstractScalarField< DGLinear<Grid> > dg ( lin.interpolate< DGLinear<Grid> >() ) ;

	ASSERT_EQ( g.nCells() * Voxel::NV, dg.size() ) ;
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
