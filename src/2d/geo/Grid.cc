/*
 * This file is part of Sand6, a C++ continuum-based granular simulator.
 *
 * Copyright 2016 Gilles Daviet <gilles.daviet@inria.fr> (Inria - Universit√© Grenoble Alpes)
 *
 * Sand6 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * Sand6 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with Sand6.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Grid.hh"
#include "BoundaryInfo.hh"

namespace d6 {

GridIterator& GridIterator::operator ++()
{
	++cell[1] ;
	if(cell[1] == grid.dim()[1]) {
		cell[1] = 0 ;
		++cell[0] ;
	}

	return *this ;
}

Index GridIterator::index() const
{
	return grid.cellIndex( cell ) ;
}

Grid::Grid(const Vec &box, const VecWi &res, const Particles *)
    : Base()
{
	m_dim = res ;
	set_box( box ) ;
}

void Grid::set_box( const Vec& box )
{
	m_dx.array() = box.array()/m_dim.array().cast< Scalar >() ;
}

Index Grid::nEdges() const
{
	return
	        m_dim[0] * (m_dim[1]+1)  +
	        m_dim[1] * (m_dim[0]+1)  ;
}

Index Grid::nFaces() const
{
	return nEdges() ;
}

void Grid::boundaryInfo( const Location &loc, const BoundaryMapper& mapper, BoundaryInfo &info ) const
{
	constexpr Scalar eps = 1.e-6 ;

	const Vec &p = pos( loc ) ;
	const Vec &b = box() ;

	info.bc = BoundaryInfo::Interior ;

	if( p[0] < eps )
		info.combine( mapper( "left"), Vec(-1,0) ) ;
	if( p[0] > b[0] - eps )
		info.combine( mapper("right"), Vec( 1,0) ) ;

	if( p[1] < eps )
		info.combine( mapper("bottom"), Vec(0,-1) ) ;
	if( p[1] > b[1] - eps )
		info.combine( mapper(   "top"), Vec(0, 1) ) ;
}

void Grid::each_neighbour( const Cell& cell, const std::function<void(const Cell&) >& func ) const
{

	Grid::Cell min = cell - Grid::Cell::Ones() ;
	Grid::Cell max = cell + Grid::Cell::Ones() ;
	clamp_cell( min );
	clamp_cell( max );

	Grid::Cell nb ;

	for( nb[0] = min[0] ; nb[0] <= max[0] ; ++nb[0] ) {
		for( nb[1] = min[1] ; nb[1] <= max[1] ; ++nb[1] ) {
			func(nb) ;
		}
	}
}

} //d6
