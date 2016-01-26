
#include "BoundaryInfo.hh"
#include "Tensor.hh"

#include "string.hh"

#include <Eigen/Geometry>

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
	case BoundaryInfo::Corner:
		proj = normal.normalized() * normal.normalized().transpose() ;
		break ;
	case BoundaryInfo::Interior:
	case BoundaryInfo::Free:
		proj.setIdentity() ;
		break ;
	}
}

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

void BoundaryInfo::combine( const Bc bc_, const Vec n )
{
	*this = combine( *this, BoundaryInfo( bc_, n ) ) ;
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

void BoundaryInfo::spinProj(MatR &proj) const
{
	proj.setIdentity() ;
}

StrBoundaryMapper::StrBoundaryMapper(const std::string &str)
{
	std::string parse_str ;

	//Predefs
	if( str == "cuve") {
		parse_str = "top:normal left:slip right:slip front:slip back:slip bottom:stick" ;
	} else if( str == "poiseuille") {
		parse_str = "top:normal left:stick right:tick front:slip back:slip bottom:normal" ;
	} else
		parse_str = str ;


	std::istringstream in( parse_str ) ;
	std::string line ;
	std::vector< std::string > tok ;

	while( in >> line ) {
		tok.clear() ;
		split( line, ":", tok );
		if( tok.size() == 2 ) {
			m_bc[canonicalize(tok[0])] = from_string(canonicalize(tok[1])) ;
		}
	}
}

BoundaryInfo::Bc StrBoundaryMapper::operator ()( const std::string &domain ) const
{
	Map::const_iterator it = m_bc.find( domain ) ;
	if( it == m_bc.end() )
		return BoundaryInfo::Stick ;

	return it->second ;
}

BoundaryInfo::Bc StrBoundaryMapper::from_string(const std::string &bc)
{
	if( bc == "slip"   ) return BoundaryInfo::Slip   ;
	if( bc == "free"   ) return BoundaryInfo::Free   ;
	if( bc == "normal" ) return BoundaryInfo::Normal ;
	return BoundaryInfo::Stick ;
}


} //ns hyb2d
