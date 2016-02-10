
#include "utils/Config.hh"

#include "geo/MeshShapeFunction.hh"
#include "geo/P2ShapeFunction.hh"
#include "geo/TetGrid.hh"

#include "ScalarField.hh"

#include <gtest/gtest.h>

using namespace d6 ;

TEST( geo, p2 )
{
	VecWi dim(2,2) ;
	Vec   box(2.,2.) ;

	TetGrid g ( box,   dim ) ;
	TetGrid gd( box, 2*dim ) ;

	ASSERT_EQ(  8, g.nCells() ) ;
	ASSERT_EQ( 16, g.nEdges() ) ;

	P2< TetGrid > p2( g ) ;
	ASSERT_EQ( 25, p2.nDOF() ) ;
	ASSERT_EQ( gd.nNodes(), p2.nDOF() ) ;

	AbstractScalarField< P2< TetGrid > > f( p2 ) ;
	f.set_constant( 1. ) ;
	ASSERT_DOUBLE_EQ( 1., f(Vec(.3,.3)) ) ;
	ASSERT_DOUBLE_EQ( 1., f(Vec(1.3,.3)) ) ;
	ASSERT_DOUBLE_EQ( 1., f(Vec(.3,1.3)) ) ;
	ASSERT_DOUBLE_EQ( 1., f(Vec(1.3,1.3)) ) ;


	f.add_at( g.locate(Vec(.5,.5)), 1. );
	ASSERT_DOUBLE_EQ( 2., f(Vec(.5,.5)) ) ;
	ASSERT_DOUBLE_EQ( 1., f(Vec(1.,1.)) ) ;

	for( auto qpIt = p2.qpBegin() ; qpIt != p2.qpEnd() ; ++qpIt ) {
		TetGrid::Location loc, loc_bis ;
		qpIt.locate(loc) ;
		p2.locate( qpIt.pos(), loc_bis ) ;

		ASSERT_TRUE( loc.cell.isApprox( loc_bis.cell ) ) ;
		ASSERT_TRUE( loc.coords.isApprox( loc_bis.coords ) ) ;

	}

}
