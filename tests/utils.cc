
#include "utils/alg.hh"

#include <gtest/gtest.h>

using namespace d6 ;

TEST( utils, alg )
{
	constexpr Index D = 6 ;
	const Index m = 5 ;

	DynVec s(m) ;
	for( int i = 0 ; i < m ; ++i ) s[i] = i+1 ;

	DynVec t ;
	set_compwise<D>( t, s ) ;

	for( int i = 0 ; i < m ; ++i ) {
		EXPECT_TRUE( Segmenter<D>::segment(t,i).isApprox( Vec6::Constant(s[i]) ) ) ;
	}

	component<D>(t, 0) *= 2 ;

	for( int i = 0 ; i < m ; ++i ) {
		EXPECT_DOUBLE_EQ( 2.*s[i], Segmenter<D>::segment(t,i)[0] ) ;
	}

	s = component<D>((const DynVec&) t, 0 ) ;
	div_compwise<D>( t, s ) ;

	for( int i = 0 ; i < m ; ++i ) {
		EXPECT_DOUBLE_EQ( 1., Segmenter<D>::segment(t,i)[0] ) ;
		EXPECT_DOUBLE_EQ( (D+1.)/2, Segmenter<D>::segment(t,i).sum() ) ;
	}

	mul_compwise<D>( t, s ) ;
	for( int i = 0 ; i < m ; ++i ) {
		EXPECT_DOUBLE_EQ( s[i], Segmenter<D>::segment(t,i)[0] ) ;
	}
}
