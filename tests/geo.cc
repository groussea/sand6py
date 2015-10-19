
#include <gtest/gtest.h>

#include "geo/Grid.hh"

#include "geo/ScalarField.hh"

using namespace d6 ;

TEST( geo, grid )
{
	Vec3i dim( 10, 5, 4 ) ;
	Vec   box( 1, 1, 1 ) ;
	Grid g( box, dim ) ;

	ASSERT_EQ( 200, g.nCells() ) ;
	ASSERT_EQ( 11*6*5, g.nNodes() ) ;
	ASSERT_TRUE( box.isApprox( g.box() ) ) ;

	Grid::Location loc ;
	g.locate( Vec(0,0,0), loc );

	ASSERT_EQ( 0, loc.cidx ) ;
	ASSERT_EQ( 0, loc.nodes[0] ) ;
	ASSERT_EQ( dim[2]+1, loc.nodes[2] ) ;
	ASSERT_EQ( (dim[1]+1)*(dim[2]+1), loc.nodes[4] ) ;
	ASSERT_EQ( (dim[1]+1)*(dim[2]+1)+(dim[2]+1)+1, loc.nodes[7] ) ;
	ASSERT_DOUBLE_EQ( 1, loc.coeffs[0] ) ;

	g.locate( Vec(1,1,1), loc );

	ASSERT_EQ( g.nCells()-1, loc.cidx ) ;
	ASSERT_EQ( g.nNodes()-1, loc.nodes[7] ) ;
	ASSERT_DOUBLE_EQ( 1, loc.coeffs[7] ) ;
}

TEST( geo, field )
{
	Vec3i dim( 10, 10, 10 ) ;
	Vec   box( 1, 1, 1 ) ;
	Grid g( box, dim ) ;

	ScalarField< Grid > phi( g ) ;
	ASSERT_EQ( g.nNodes(), phi.flatten().rows() ) ;

}

