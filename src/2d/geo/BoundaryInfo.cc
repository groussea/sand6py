
#include "BoundaryInfo.hh"
#include "Tensor.hh"

#include "string.hh"

#include <Eigen/Geometry>

namespace d6 {

BoundaryInfo BoundaryInfo::combine( const BoundaryInfo& b1, const BoundaryInfo& b2 )
{
	if( b1.bc == Interior )	return b2 ;
	if( b2.bc == Interior )	return b1 ;
	if( b1.bc == Free )		return b2 ;
	if( b2.bc == Free )		return b1 ;

	if( b1.bc == Stick )	return b1 ;
	if( b2.bc == Stick )	return b2 ;

	const Scalar nd = std::fabs( b1.normal.dot(b2.normal) ) ;

	if( b1.bc == b2.bc ||
			(b1.bc == Normal && b2.bc == Corner) ||
			(b1.bc == Corner && b2.bc == Normal)
			) {
		// Same constraints
		if( std::fabs( 1-nd ) < 1.e-6 )
			return b1 ;

		// Incompatible free direction
		return BoundaryInfo( Stick, b1.normal ) ;
	} else {

		// Incompatible free direction
		if( nd > 1.e-6 )
			return BoundaryInfo( Stick, b1.normal ) ;

		if( b1.bc == Normal || b1.bc == Corner )
			return b1 ;

		return b2 ;
	}

}


void BoundaryInfo::stressProj( MatS &proj ) const
{
	proj.setIdentity() ;

	switch( bc )
	{
	case BoundaryInfo::Corner:
	case BoundaryInfo::Interior:
	case BoundaryInfo::Stick:
	case BoundaryInfo::Free:
		return ;
	default:
		break ;
	}

	const Vec & n = normal ;

	Eigen::Matrix< Scalar, WD, SD > N ;
		// [ n0 ;  n0 ; n1 ] . taubar = 0
		// [ n1 ; -n1 ; n0 ] . taubar = 0
	N << n[0],  n[0], n[1],
		 n[1], -n[1], n[0];

	N.row(0) = N.row(0).normalized() ;
	N.row(1) = N.row(1).normalized() ;

	switch( bc )
	{
	case BoundaryInfo::Free:
		// (\tau n) = 0
		break ;
	case BoundaryInfo::Normal:
		// nn' (\tau n) = 0
		N = (n*n.transpose()) * N ;
		break ;
	case BoundaryInfo::Slip:
		// (\tau n) - nn' (\tau n) = 0
		N = N - (n*n.transpose()) * N ;
		break ;
	default:
		break ;
	}

	proj -=	  N.row(0).transpose()*N.row(0)
			+ N.row(1).transpose()*N.row(1) ;
}

} //ns hyb2d
