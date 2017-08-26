#include "MeshShapeFunction.hh"

#include "Grid.hh"

namespace d6 {

template<>
void Facewise<Grid>::list_nodes( const Location& loc, typename Base::NodeList& nodes ) const
{
	const Grid& grid = Base::mesh().derived() ;

	int off = 0 ;
	for( int j = 0 ; j < 3 ; ++j )
	{
		const int k = (j+1)%3 ;
		const int l = (j+2)%3 ;
		const Index nk = grid.dim()[k] ;
		const Index nl = grid.dim()[l] ;

		nodes[2*j+0] = off + loc.cell[j]*nk*nl + loc.cell[k]*nl + loc.cell[l] ;
		nodes[2*j+1] = off + (loc.cell[j]+1)*nk*nl + loc.cell[k]*nl + loc.cell[l] ;
		off += grid.dim()[j] * nk * nl ;
	}
}


template<>
void Facewise<Grid>::dof_coeffs( const typename MeshType::Coords& coords, typename Base::CoefList& coeffs ) const
{
	for( int j = 0 ; j < 3 ; ++j )
	{
		coeffs[2*j + 0] = (1-coords[j]) ;
		coeffs[2*j + 1] = (  coords[j]) ;
	}
}

template<>
void Facewise<Grid>::get_derivatives( const Location&, typename Base::Derivatives& dc_dx ) const
{
	dc_dx.setZero() ;

	for( int j = 0 ; j < 3 ; ++j )
	{
		dc_dx(2*j+0, j) = -1/mesh().dx()[j] ;
		dc_dx(2*j+1, j) =  1/mesh().dx()[j] ;
	}
}

template<>
void Edgewise<Grid>::list_nodes( const Location& loc, typename Base::NodeList& nodes ) const
{
	const Grid& grid = Base::mesh().derived() ;

	int off = 0 ;
	for( int j = 0 ; j < 3 ; ++j )
	{
		const int k = (j+1)%3 ;
		const int l = (j+2)%3 ;
		const Index nk = grid.dim()[k] + 1;
		const Index nl = grid.dim()[l] + 1;

		nodes[4*j+0] = off + loc.cell[j]*nk*nl + loc.cell[k]*nl + loc.cell[l] ;
		nodes[4*j+1] = off + loc.cell[j]*nk*nl + loc.cell[k]*nl + loc.cell[l]+1 ;
		nodes[4*j+2] = off + loc.cell[j]*nk*nl + (loc.cell[k]+1)*nl + loc.cell[l] ;
		nodes[4*j+3] = off + loc.cell[j]*nk*nl + (loc.cell[k]+1)*nl + loc.cell[l]+1 ;
		off += grid.dim()[j] * nk * nl ;
	}
}


template<>
void Edgewise<Grid>::dof_coeffs( const typename MeshType::Coords& coords, typename Base::CoefList& coeffs ) const
{
	for( int j = 0 ; j < 3 ; ++j )
	{
		const int k = (j+1)%3 ;
		const int l = (j+2)%3 ;

		coeffs[4*j + 0] = (1-coords[k])*(1-coords[l]) ;
		coeffs[4*j + 1] = (1-coords[k])*(coords[l]) ;
		coeffs[4*j + 2] = (coords[k])  *(1-coords[l]) ;
		coeffs[4*j + 3] = (coords[k])  *(coords[l]) ;
	}
}

template<>
void Edgewise<Grid>::get_derivatives( const Location& loc, typename Base::Derivatives& dc_dx ) const
{
	const typename MeshType::Coords& coords = loc.coords ;

	for( int j = 0 ; j < 3 ; ++j )
	{
		const int k = (j+1)%3 ;
		const int l = (j+2)%3 ;

		dc_dx(4*j+0, j) = 0 ;
		dc_dx(4*j+0, k) = -(1-coords[l]) ;
		dc_dx(4*j+0, l) = -(1-coords[k]) ;

		dc_dx(4*j+1, j) = 0 ;
		dc_dx(4*j+1, k) = -(coords[l]) ;
		dc_dx(4*j+1, l) = (1-coords[k]) ;

		dc_dx(4*j+2, j) = 0 ;
		dc_dx(4*j+2, k) = (1-coords[l]) ;
		dc_dx(4*j+2, l) = -(coords[k]) ;

		dc_dx(4*j+3, j) = 0 ;
		dc_dx(4*j+3, k) = (coords[l]) ;
		dc_dx(4*j+3, l) = (coords[k]) ;
	}

	for (int k = 0 ; k < WD ; ++k)
		dc_dx.col( k ) /= mesh().dx()[k] ;
}

} // d6
