
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

void BoundaryInfo::combine( const Bc bc_, const Vec n )
{
	*this = combine( *this, BoundaryInfo( bc_, n ) ) ;
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
