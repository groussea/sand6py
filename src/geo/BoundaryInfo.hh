#ifndef D6_BOUNDARY_INFO_HH
#define D6_BOUNDARY_INFO_HH

#include "utils/alg.hh"

#include <string>
#include <vector>

namespace d6 {

struct BoundaryInfo {

	enum Bc {
		Interior,
		Stick,   // u = 0
		Slip,    // u.n = 0, d( (I - nn') u ).n = 0, (I - nn')sn = 0
		Normal,  // (I -nn') u = 0, d ( u.n ).n = 0, (nn')sn = 0
		Free     // d( u ).n = 0, sn = 0
	};

	Bc bc ;
	Vec normal ;

	BoundaryInfo() :
		bc( Interior )
	{}

	void set( const Bc bc_, const Vec n )
	{
		bc = bc_ ;
		normal = n ;
	}

	void    velProj( Mat &proj ) const ;
	void   spinProj( Mat &proj ) const ;
	void stressProj( Mat66 &proj ) const ;

};

typedef std::vector< BoundaryInfo > BoundaryConditions ;

struct BoundaryMapper {
	virtual BoundaryInfo::Bc operator() ( const std::string &/*domain*/ ) const
	{ return BoundaryInfo::Stick ; }
};

} //ns d6

#endif
