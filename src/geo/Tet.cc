#include "Tet.hh"

#include "Tensor.hh"

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

void Tet::update_geometry( unsigned char symmetry, int num )
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

	geometry.sym ^= symmetry ;

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

void Tet::local_coords(const Vec &pos, Coords &coords) const
{
	const Vec scaled = pos.array()/box ;

	coords[0] = 1 - scaled[2] - scaled[geometry.part] ;
	coords[1] = scaled[geometry.part] - scaled[1-geometry.part] ;
	coords[2] = scaled[1-geometry.part] ;
	coords[3] = scaled[2] ;
}

void Tet::compute_derivatives( const Coords & coords, Derivatives& dc_dx ) const
{
	(void) coords ;

	dc_dx.setZero() ;
	dc_dx(0,  geometry.part) = -1./box[  geometry.part] ;
	dc_dx(1,  geometry.part) =  1./box[  geometry.part] ;
	dc_dx(1,1-geometry.part) = -1./box[1-geometry.part] ;
	dc_dx(2,1-geometry.part) =  1./box[1-geometry.part] ;
	dc_dx(0,2) = -1/box[2] ;
	dc_dx(3,2) =  1/box[2] ;

	//to_word

	if( geometry.rotate(0) ) {
		dc_dx.col(2).swap( dc_dx.col(1) ) ;
		dc_dx.col(1) *= -1 ;
	}
	if( geometry.rotate(1) ) {
		dc_dx.col(2).swap( dc_dx.col(0) ) ;
		dc_dx.col(2) *= -1 ;
	}
	if( geometry.rotate(2) ) {
		dc_dx.col(1).swap( dc_dx.col(0) ) ;
		dc_dx.col(0) *= -1 ;
	}

	dc_dx.col(0) = geometry.sign(0)*dc_dx.col(0) ;
	dc_dx.col(1) = geometry.sign(1)*dc_dx.col(1) ;
	dc_dx.col(2) = geometry.sign(2)*dc_dx.col(2) ;
}

Index Tet::sample_uniform( const unsigned N, const Index start, Points &points, Frames &frames ) const
{
	/*
	const Scalar a = 1./6 ;
	const Scalar b = 1./2 ;

	(void) N ;
	const Vec subBox = box.array() / std::pow(24, 1./3) ; //Nsub.array().cast< Scalar >() ;

	Vec6 frame ;
	tensor_view( frame ).set_diag( Vec( .25 * subBox.array() * subBox.array() ) ) ;

	Index p = start ;
	for( Index k = 0 ; k < 4 ; ++k ) {
		Coords coords = Coords::Constant(a) ;
		coords[k] = b ;
		points.col(p) = pos(coords) ;
		frames.col(p) = frame ;
		++p ;
	}
	*/

	Scalar min = box.minCoeff() ;

	Vec3i Nsub ;
	for( int k = 0 ; k < 3 ; ++ k)
		Nsub[k] = N * std::round( box[k] / min ) ;

	const Vec subBox = box.array() / Nsub.array().cast< Scalar >() ;

	Vec rot_box = subBox ;
	to_world( rot_box ) ;

	Vec6 frame ;
	tensor_view( frame ).set_diag( Vec( .25 * rot_box.array() * rot_box.array() ) ) ;

	Index p = start ;

	Vec local ;
	Coords coords ;

	for( int i = 0 ; i < Nsub[0] ; ++i )
		for( int j = 0 ; j < Nsub[1] ; ++j )
			for( int k = 0 ; k < Nsub[2] ; ++k ) {

				local = (Vec(i+.5,j+.5,k+.5).array() * subBox.array()).matrix() ;
				local_coords( local, coords );
				if( coords.minCoeff() < 0 )
					continue ;

				points.col(p) = pos( coords )  ;
				frames.col(p) = frame ;
				++p ;
			}

	return p - start;
}

} //d6
