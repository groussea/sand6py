#include "Tet.hh"

namespace d6 {


void Tet::to_world( Vec& pos ) const {

	if( orientation & 2 ) {
		pos = Vec( pos[0], -pos[2], pos[1] ) ;
	}
	if( orientation & 4 ) {
		pos = Vec( pos[2], pos[1], -pos[0] ) ;
	}
	if( orientation & 8 ) {
		pos = Vec( -pos[1], pos[0], pos[2] ) ;
	}
}
void Tet::to_local( Vec& pos ) const {
	if( orientation & 8 ) {
		pos = Vec( pos[1], -pos[0], pos[2] ) ;
	}
	if( orientation & 4 ) {
		pos = Vec( -pos[2], pos[1], pos[0] ) ;
	}
	if( orientation & 2 ) {
		pos = Vec( pos[0], pos[2], -pos[1] ) ;
	}
}
void Tet::offset( int cornerIndex, Vec &v ) const
{

	v.setZero() ;

	switch (cornerIndex) {
	case 1:
		v[orientation&1] += box[orientation&1] ;
		break;
	case 2:
		v[0] = box[0] ;
		v[1] = box[1] ;
		break ;
	case 3:
		v[2] = box[2] ;
		break ;
	default:
		break ;
	}
}

void Tet::compute_vertices( Vertices& vertices ) const
{
	vertices.setZero() ;
	vertices(orientation&1, 1) = box[orientation&1] ;
	vertices(0,2) = box[0] ;
	vertices(1,2) = box[1] ;
	vertices(2,3) = box[2] ;
}

Tet::QuadPoints Tet::Qps()
{

	const Scalar a = (5. -   std::sqrt(5.) ) / 20 ;
	const Scalar b = (5. + 3*std::sqrt(5.) ) / 20 ;

	QuadPoints qps ;
	qps.setConstant( a ) ;
	qps.diagonal().setConstant( b ) ;

	return qps ;
}

void Tet::local_coords(const Vec &pos, Coords &coords) const
{
	const Vec scaled = pos.array()/box ;
	const int dir = orientation&1 ;


	coords[3] = scaled[2] ;
	coords[2] = scaled[1-dir] ;
	coords[1] = scaled[dir] - coords[2] ;
	coords[0] = 1 - coords.segment<3>(1).sum() ;
}

} //d6
