#include "Tet.hh"

namespace d6 {

void Tet::to_world( Vec& pos ) const {

	if( geometry.rotate(0) ) {
		pos = Vec( pos[0], -pos[2], pos[1] ) ;
	}
	if( geometry.rotate(1) ) {
		pos = Vec( pos[2], pos[1], -pos[0] ) ;
	}
	if( geometry.rotate(2) ) {
		pos = Vec( -pos[1], pos[0], pos[2] ) ;
	}

	pos[0] = geometry.sign(0)*pos[0] ;
	pos[1] = geometry.sign(1)*pos[1] ;
	pos[2] = geometry.sign(2)*pos[2] ;

}
void Tet::to_local( Vec& pos ) const {

	pos[0] = geometry.sign(0)*pos[0] ;
	pos[1] = geometry.sign(1)*pos[1] ;
	pos[2] = geometry.sign(2)*pos[2] ;

	if( geometry.rotate(2) ) {
		pos = Vec( pos[1], -pos[0], pos[2] ) ;
	}
	if( geometry.rotate(1) ) {
		pos = Vec( -pos[2], pos[1], pos[0] ) ;
	}
	if( geometry.rotate(0) ) {
		pos = Vec( pos[0], pos[2], -pos[1] ) ;
	}
}

void Tet::update_geometry( unsigned char rotation, int num )
{
	geometry.part = num %2 ;

	switch( num/2 ) {
		case 1:
			// Rot y, sym x
			geometry.rot = 1 << 1 ;
			geometry.sym = 1 << 0 ;
			break ;
		case 2:
			// Rot x, sym z
			geometry.rot = 1 << 0 ;
			geometry.sym = 1 << 2 ;
			break ;
		default:
			geometry.rot = 0 ;
			geometry.sym = 0 ;
			break ;
	}

	geometry.rot ^= rotation ;

	Vec o (-.5,-.5,-.5)	 ;
	Vec o_rot = o ;
	to_world( o_rot );

	origin.array() += box*( o_rot - o ).array() ;

	Vec rot_box = box ;
	to_local( rot_box ) ;
	box = rot_box.array().abs() ;

}



void Tet::offset( int cornerIndex, Vec &v ) const
{

	v.setZero() ;

	switch (cornerIndex) {
		case 1:
			v[geometry.part] += box[geometry.part] ;
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
	vertices(geometry.part, 1) = box[geometry.part] ;
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
	const int dir = geometry.part ;


	coords[3] = scaled[2] ;
	coords[2] = scaled[1-dir] ;
	coords[1] = scaled[dir] - coords[2] ;
	coords[0] = 1 - coords.segment<3>(1).sum() ;
}

} //d6
