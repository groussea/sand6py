#include "UnstructuredShapeFunction.hh"

#include "Grid.hh"

namespace d6 {


void UnstructuredDOFs::compute_weights_from_vertices( const Vec &box, const VecWi& res )
{

	//FIXME -- chains of particles -> max density

	const Index n = count() ;
	weights.setConstant( vertices.cols(), 0 ) ;

	DynArr dist = DynArr::Ones( n ) ;

	Grid g( box, res ) ;
	std::vector< std::vector< Index > > ids( g.nCells() ) ;

	for( Index i = 0 ; i < n ; ++i ) {
		Grid::Location loc ;
		g.locate( vertices.col(i), loc) ;
		ids[ g.cellIndex(loc.cell) ].push_back( i ) ;
	}

#pragma omp parallel for
	for( Index i = 0 ; i < n ; ++i ) {
		Grid::Location loc ;
		g.locate( vertices.col(i), loc) ;

		for( Index k = -1 ; k <2  ; ++k ) {
			for( Index j = -1 ; j< 2 ; ++j ) {
				Grid::Cell nb = loc.cell ;
				nb[0] += j ; nb[1] += k ;
				g.clamp_cell( nb ); //TODO optimize

				for ( Index pid : ids[g.cellIndex(nb)] ) {
					if( pid == i ) continue ;
					Scalar d2 = ( vertices.col(i) - vertices.col(pid) ).squaredNorm() ;
					if( d2 < dist(i) )
						dist(i) = d2 ;
				}
			}
		}
	}

	weights = dist ; // Assume sqaure shape

}

} //d6
