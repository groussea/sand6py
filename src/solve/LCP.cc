#include "LCP.hh"

#include "utils/Log.hh"

#include <Eigen/Eigenvalues>
#include <bogus/Core/Block.impl.hpp>
#include <bogus/Core/BlockSolvers.impl.hpp>

namespace d6 {

LCP::LCP(const LCPData &data)
	: m_data( data )
{}

static void ackResidual( unsigned iter, Scalar res ) {
	Log::Debug() << "LCP " << iter << " =\t " << res << std::endl ;
}


Scalar LCP::solve(DynVec &x) const
{
	typedef typename bogus::SparseBlockMatrix< Scalar, bogus::SYMMETRIC > WType ;

	WType W = m_data.H * m_data.H.transpose() ;

	bogus::GaussSeidel< WType > gs( W ) ;
	gs.callback().connect( ackResidual );

	return gs.solve( bogus::LCPLaw< Scalar>(), m_data.w, x ) ;
}

} //d6
