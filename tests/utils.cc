
#include <gtest/gtest.h>

#include "utils/Config.hh"

using namespace d6 ;

TEST( utils, config )
{
	Config c ;
	ASSERT_TRUE( c.from_file("../scenes/test.conf") ) ;

	ASSERT_EQ( c.res, Vec3i(1, 5, 10) ) ;
	ASSERT_TRUE( c.gravity.isApprox( Vec(40,0,0) ) ) ;
	ASSERT_TRUE( c.box.isApprox( Vec(10,20,30) ) ) ;
	ASSERT_DOUBLE_EQ( c.dt(), 1.e-2 ) ;

	const Scalar Re = ( c.box.minCoeff() * std::sqrt( c.box.minCoeff() * c.gravity.norm()) * c.volMass ) / c.viscosity ;
	ASSERT_DOUBLE_EQ( 2.e8, Re ) ;

	c.internalize();

	ASSERT_TRUE( c.box.isApprox( Vec(1,2,3) ) ) ;
	ASSERT_TRUE( c.gravity.isApprox( Vec(1,0,0) ) ) ;
	ASSERT_DOUBLE_EQ( c.units().U, 20. ) ;
	ASSERT_DOUBLE_EQ( c.units().T, .5 ) ;
	ASSERT_DOUBLE_EQ( c.units().P, 4.e5 ) ;
	ASSERT_DOUBLE_EQ( c.units().M, 2.e5 ) ;
	ASSERT_DOUBLE_EQ( c.units().M, c.units().R * c.units().U * c.units().L ) ;

	ASSERT_DOUBLE_EQ( c.viscosity, 1./Re ) ;

}
