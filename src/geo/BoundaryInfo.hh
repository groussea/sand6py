#ifndef D6_BOUNDARY_INFO_HH
#define D6_BOUNDARY_INFO_HH

#include "utils/alg.hh"

#include <string>
#include <vector>

#include <unordered_map>

namespace d6 {

struct BoundaryInfo {

	enum Bc {
		Interior,
		Stick,   // u = 0
		Slip,    // u.n = 0, d( (I - nn') u ).n = 0, (I - nn')sn = 0
		Normal,  // (I -nn') u = 0, d ( u.n ).n = 0, (nn')sn = 0
		Free,    // d( u ).n = 0, sn = 0
		Corner
	};

	Bc bc ;
	Vec normal ;

	BoundaryInfo() :
		bc( Interior )
	{}

	BoundaryInfo( const Bc bc_, const Vec n )
		: bc(bc_), normal(n)
	{}

	void set( const Bc bc_, const Vec n )
	{
		bc = bc_ ;
		normal = n ;
	}

	void combine( const Bc bc_, const Vec n ) ;
	static BoundaryInfo combine(const BoundaryInfo &b1, const BoundaryInfo &b2 ) ;

	void    velProj( Mat  &proj ) const ;
	void   spinProj( MatR &proj ) const ;
	void stressProj( MatS &proj ) const ;

};

typedef std::vector< BoundaryInfo > BoundaryConditions ;

struct BoundaryMapper {
	virtual BoundaryInfo::Bc operator() ( const std::string &/*domain*/ ) const
	{ return BoundaryInfo::Stick ; }
};

class StrBoundaryMapper : public BoundaryMapper
{
public:

	explicit StrBoundaryMapper( const std::string & str ) ;

	virtual BoundaryInfo::Bc operator() ( const std::string &domain ) const ;

private:

	BoundaryInfo::Bc from_string( const std::string &bc ) ;

	typedef std::unordered_map< std::string, BoundaryInfo::Bc > Map ;
	Map m_bc ;
};

} //ns d6

#endif
