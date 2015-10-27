
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

	switch( bc )
	{
	case BoundaryInfo::Free:
	case BoundaryInfo::Normal:
		// (\tau n n') n = \tau n (n'n) = \tau n = (\tau I) n )
		// (\tau I ) n = 0
		// [ sqrt2_3 n0 ;  n0 ; -1/sqrt3 n0 ; n1 ; n2 ;  0 ] . taubar = 0
		// [ sqrt2_3 n1 ; -n1 ; -1/sqrt3 n1 ; n0 ;  0 ; n2 ] . taubar = 0
		// [ sqrt2_3 n2 ;   0 ;  2/sqrt3 n2 ;  0 ; n0 ; n1 ] . taubar = 0
	{
		const Vec & n = normal ;
		Vec6 N1 ; N1 << s_sqrt23 * n[0],  n[0],  -s_isqrt3 * n[0], n[1], n[2],   0  ;
		Vec6 N2 ; N2 << s_sqrt23 * n[1], -n[1],  -s_isqrt3 * n[1], n[0],   0 , n[2] ;
		Vec6 N3 ; N3 << s_sqrt23 * n[2],    0 , 2*s_isqrt3 * n[2],   0 , n[0], n[1] ;

		proj -= N1*N1.transpose() + N2*N2.transpose() + N3*N3.transpose() ;
	}
	case BoundaryInfo::Interior:
	case BoundaryInfo::Stick:
	case BoundaryInfo::Slip:
		//		// (\tau ( I - n n' ) ) n = \tau n - \tau n (n'n) = 0
		break ;
	}
}

void BoundaryInfo::spinProj(Mat &proj) const
{
	proj.setIdentity() ;
}


} //ns hyb2d
