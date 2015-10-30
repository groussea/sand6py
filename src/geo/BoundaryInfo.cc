
#include "BoundaryInfo.hh"

namespace d6 {

void BoundaryInfo::velProj( Mat &proj ) const
{
	switch( bc )
	{
	case BoundaryInfo::Stick:
		proj.setZero() ;
		break ;
	case BoundaryInfo::Slip:
		proj = Mat::Identity() - normal * normal.transpose() ;
		break ;
	case BoundaryInfo::Normal:
		proj = normal * normal.transpose() ;
		break ;
	case BoundaryInfo::Interior:
	case BoundaryInfo::Free:
		proj.setIdentity() ;
		break ;
	}
}

void BoundaryInfo::stressProj( Mat66 &proj ) const
{
	static constexpr Scalar s_sqrt23 = std::sqrt(2.)/std::sqrt(3.) ;
	static constexpr Scalar s_isqrt3 = 1./std::sqrt(3.) ;

	proj.setIdentity() ;

	/*
	 * nn' tau =
	 * 00 : ( s23 a + b - c3 ) n0^2 + d n0n1 + e n0n2
	 * 01 : ( s23 a - b - c3 ) n0n1 + d n0^2 + f n0n2
	 * 02 : ( s23 a     +2c3 ) n0n2 + e n0^2 + f n0n1
	 * 10 : ( s23 a + b - c3 ) n0n1 + d n1^2 + e n1n2
	 * 11 : ( s23 a - b - c3 ) n1^2 + d n0n1 + f n1n2
	 * 12 : ( s23 a     +2c3 ) n1n2 + e n0n1 + f n1^2
	 * 20 : ( s23 a + b - c3 ) n0n2 + d n1n2 + e n2^2
	 * 21 : ( s23 a - b - c3 ) n1n2 + d n0n2 + f n2^2
	 * 22 : ( s23 a     +2c3 ) n2^2 + e n0n2 + f n1n2
	 *
	 * nn' tau n
	 * 0 : n0( ( s23 a + b - c3 ) n0^2 + d n0n1 + e n0n2 ) + n1 (( s23 a - b - c3 ) n0n1 + d n0^2 + f n0n2 ) + n2(( s23 a     +2c3 ) n0n2 + e n0^2 + f n0n1)
	 *   = s23 a ( n0^3 + n0n1^2 + n0 n2^2 ) + b ( n0^3 - n0n1^2 ) + c3 ( 2 n0n2^2 - n0^3 - n0n1^2 ) +d ( 2 n1n0^2 ) + e ( 2 n0^2n2 ) + f( 2n0n1n2 )
	 *   = s23 a n0 + b n0 ( n0^2 - n1^2 ) + c3 n0 ( 3n2^2 - 1 ) +d ( 2 n1n0^2 ) + e ( 2 n0^2n2 ) + f( 2n0n1n2 )
	 *   /n0 = s23 a + b( n0^2 - n1^2 ) + c3 ( 3n2^2 - 1 ) + 2 d n0n1 + 2 e n0n2 + 2 f n1n2
	 *
	 * 1 : n0( ( s23 a + b - c3 ) n0n1 + d n1^2 + e n1n2 ) + n1 ( ( s23 a - b - c3 ) n1^2 + d n0n1 + f n1n2 ) + n2 ( ( s23 a     +2c3 ) n1n2 + e n0n1 + f n1^2 )
	 *   /n1 = s23 a( n0^2 + n1^2 + n2^2 ) + b ( n0^2 - n1^2 ) + c3 ( 2n2^2 -n0^2 -n1^2 ) +d (2n0n1 ) + e(2 n2n0) + f (2 n1 n2)
	 *       = s23 a + b( n0^2 - n1^2 ) + c3 ( 3n2^2 - 1 ) + 2 d n0n1 + 2 e n0n2 + 2 f n1n2
	 *
	 *  2: ....
	 */

	const Vec & n = normal ;
	Vec6 N1, N2, N3 ;
	switch( bc )
	{
	case BoundaryInfo::Free:
		// (\tau n) = 0
		// [ sqrt2_3 n0 ;  n0 ; -1/sqrt3 n0 ; n1 ; n2 ;  0 ] . taubar = 0
		// [ sqrt2_3 n1 ; -n1 ; -1/sqrt3 n1 ; n0 ;  0 ; n2 ] . taubar = 0
		// [ sqrt2_3 n2 ;   0 ;  2/sqrt3 n2 ;  0 ; n0 ; n1 ] . taubar = 0
	{
		N1 << s_sqrt23 * n[0],  n[0],  -s_isqrt3 * n[0], n[1], n[2],   0  ;
		N2 << s_sqrt23 * n[1], -n[1],  -s_isqrt3 * n[1], n[0],   0 , n[2] ;
		N3 << s_sqrt23 * n[2],    0 , 2*s_isqrt3 * n[2],   0 , n[0], n[1] ;

		proj -= N1*N1.transpose() + N2*N2.transpose() + N3*N3.transpose() ;
		break ;
	}
	case BoundaryInfo::Normal:
		// nn' (\tau n) = 0
		N1 << s_sqrt23,  n[0]*n[0]-n[1]*n[1], s_isqrt3 * (3*n[2]*n[2]-1), 2*n[0]*n[1], 2*n[0]*n[2],   2*n[1]*n[2]  ;
		proj -= N1*N1.transpose()  ;
		break ;
	case BoundaryInfo::Slip:
		// (\tau n) - nn' (\tau n) = 0
		N1 << 0,  n[0]*( 1 -n[0]*n[0] +n[1]*n[1]), -s_isqrt3 * n[0] * 3*n[2]*n[2], n[1] * ( 1 - 2*n[0]*n[0]), n[2] * (1 - 2*n[0]*n[0]), -2*n[0]*n[1]*n[2]  ;
		N2 << 0,  n[1]*(-1 -n[0]*n[0] +n[1]*n[1]), -s_isqrt3 * n[1] * 3*n[2]*n[2], n[0] * ( 1 - 2*n[1]*n[1]), -2*n[0]*n[1]*n[2] , n[2] * (1 - n[1]*n[1]) ;
		N3 << 0,  n[2]*(   -n[0]*n[0] +n[1]*n[1]),3*s_isqrt3 * n[2] *(1-n[2]*n[2]),-2*n[0]*n[1]*n[2] , n[0]*(1 - 2*n[2]*n[2]), n[1]*( 1 - 2*n[2]*n[2] ) ;

		// 1 : n0 ; 2 : n1 ; 3 : n2
		// n0n0 * (( 1 -n[0]*n[0] +n[1]*n[1])) + n1*n1*(-1 -n[0]*n[0] +n[1]*n[1] + n2*n2*(n[1]*n[1] -n[0]*n[0])
		// (n0^2 + n1^2 +n2^2 )*((-n[0]*n[0] +n[1]*n[1])) + n[0]*n[0] - n[1]*n[1]
		// (n0^2 + n1^2 +n2^2 -1 ) = 0

		// -( n[0]^2 + n[1]^2 ) * 3 * n[2]^2 ) + 3 *n[2]^2 * (1 - n[2]^2)
		// 3*(  -( n[0]^2 + n[1]^2 + n[2]^2  ) * n[2]^2 + n[2]^2 ) = 0

		// n[0]*( n[2] - 2n[2]n[0]^2) - 2*n[1]*n[0]n[1]n[2] + n[2]n[0]*(1 - 2*n[2]*n[2])
		// n[0]n[2]*( 1 - 2[n0]^2  - 2*n[1]n[1] + 1 - 2[n2]^2 ) = 0
		// ok rg 2

		proj -= N1*N1.transpose() + N2*N2.transpose() + N3*N3.transpose() ;
		break ;
	case BoundaryInfo::Interior:
	case BoundaryInfo::Stick:
		break ;
	}
}

void BoundaryInfo::spinProj(Mat &proj) const
{
	proj.setIdentity() ;
}


} //ns hyb2d
