
#include <gtest/gtest.h>

#include "geo/Grid.hh"

#include "geo/ScalarField.hh"
#include "geo/VectorField.hh"
#include "geo/TensorField.hh"

using namespace d6 ;

TEST( geo, grid )
{
	Vec3i dim( 10, 5, 4 ) ;
	Vec   box( 1, 1, 1 ) ;
	Grid g( box, dim ) ;

	ASSERT_EQ( 200, g.nCells() ) ;
	ASSERT_EQ( 11*6*5, g.nNodes() ) ;
	ASSERT_TRUE( box.isApprox( g.box() ) ) ;

	Grid::Interpolation itp ;
	g.interpolate( Vec(0,0,0), itp );

	ASSERT_EQ( 0, itp.nodes[0] ) ;
	ASSERT_EQ( dim[2]+1, itp.nodes[2] ) ;
	ASSERT_EQ( (dim[1]+1)*(dim[2]+1), itp.nodes[4] ) ;
	ASSERT_EQ( (dim[1]+1)*(dim[2]+1)+(dim[2]+1)+1, itp.nodes[7] ) ;
	ASSERT_DOUBLE_EQ( 1, itp.coeffs[0] ) ;

	g.interpolate( Vec(1,1,1), itp );

	ASSERT_EQ( g.nNodes()-1, itp.nodes[7] ) ;
	ASSERT_DOUBLE_EQ( 1, itp.coeffs[7] ) ;

	Index i = 0 ;
	for( GridIterator it( g.cellBegin() ) ;  it != g.cellEnd() ; ++it, ++i )
	{
		ASSERT_EQ( i, it.index() ) ;
	}
	ASSERT_EQ( i, g.nCells() ) ;
}

TEST( geo, field )
{
	Vec3i dim( 10, 10, 10 ) ;
	Vec   box( 1, 1, 1 ) ;
	Grid g( box, dim ) ;

	AbstractScalarField< Grid > phi( g ) ;
	ASSERT_EQ( g.nNodes(), phi.flatten().rows() ) ;

	phi.set_constant( 3 );
	ASSERT_DOUBLE_EQ( 3., phi( Vec( 0.2, 0.7, 0.5 ) ) ) ;

	phi.set_zero();
	ASSERT_DOUBLE_EQ( 0., phi( Vec( 0.3, 0.4, 0.1 ) ) ) ;
	phi.add_at( Vec( 0.3, 0.4, 0.1 ), 1 );
	ASSERT_FLOAT_EQ( 1,  phi( Vec( 0.3, 0.4, 0.1 ) ) ) ;
	ASSERT_FLOAT_EQ( 0.3,  phi( Vec( 0.37, 0.4, 0.1 ) ) ) ;
	ASSERT_FLOAT_EQ( 0.7,  phi( Vec( 0.3, 0.37, 0.1 ) ) ) ;
	ASSERT_FLOAT_EQ( 0.5,  phi( Vec( 0.3, 0.4, 0.15 ) ) ) ;

	AbstractVectorField< Grid > u( g ) ;
	u.set_constant( Vec(0,1,0) ) ;
	ASSERT_TRUE( Vec(0,1,0).isApprox( u( Vec( 0.2, 0.7, 0.5 ) ) ) ) ;

	u.add_at( Vec( 0.35, 0.45, 0.15 ), Vec(1,0,0) );
	ASSERT_TRUE( Vec(0.125,1,0).isApprox( u( Vec( 0.3, 0.4, 0.2 ) ) ) ) ;
}

