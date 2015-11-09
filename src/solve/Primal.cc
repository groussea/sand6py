#include "Primal.hh"

#include "utils/Log.hh"

#include <bogus/Core/Block.impl.hpp>
#include <bogus/Core/BlockSolvers.impl.hpp>
#include <bogus/Extra/SecondOrder.impl.hpp>

namespace d6 {

bool PrimalData::load(const char *file)
{
	(void) file ;
	return false ;
}

bool PrimalData::dump(const char *file) const
{
	(void) file ;
	return false ;
}


Primal::Primal(const PrimalData &data)
	: m_data( data )
{}

void ackResidual( unsigned iter, Scalar res ) {
	Log::Debug() << "GS " << iter << " =\t " << res << std::endl ;
}


Scalar Primal::solve(DynVec &lambda, DynVec &gamma) const
{
	typedef typename FormMat<6,6>::SymType WType ;
	WType W = m_data.H * m_data.H.transpose() ;

	// Add rigid bodies jacobians
	bogus::SparseBlockMatrix< Mat66, bogus::SYMMETRIC > Mi ;
	Mi.setRows(1) ;
	Mi.insertBack(0,0) ;
	Mi.finalize();

	for( unsigned i = 0 ; i < m_data.jacobians.size() ; ++i ) {

		Mi.block(0) = m_data.inv_inertia_matrices.block<6,6>(0, 6*i );

		const PrimalData::JacobianType &JM = m_data.jacobians[i] ;

		W +=  JM * Mi * JM.transpose() ;
	}


	bogus::SOCLaw< 6, Scalar, true > law( m_data.n(), m_data.mu.data() ) ;

	bogus::GaussSeidel< WType > gs( W ) ;
	gs.callback().connect( ackResidual );

	const double res = gs.solve( law, m_data.w, lambda ) ;

	gamma = W*lambda + m_data.w ;

	return res ;
}


}
