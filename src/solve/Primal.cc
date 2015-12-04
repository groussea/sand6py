#include "Primal.hh"

#include "utils/Log.hh"

#include <bogus/Core/Block.impl.hpp>
#include <bogus/Core/BlockSolvers.impl.hpp>
#include <bogus/Extra/SecondOrder.impl.hpp>

#include <bogus/Core/Utils/Timer.hpp>

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

Primal::SolverOptions::SolverOptions()
	: algorithm( GaussSeidel ),
	  maxIterations(250), maxOuterIterations( 15 ),
	  tolerance( 1.e-6 )
{}


Primal::Primal(const PrimalData &data)
	: m_data( data )
{}

void Primal::SolverStats::ackResidual( unsigned iter, Scalar res )
{
	nIterations = iter ;
	Log::Debug() << "Primal " << iter << " =\t " << res << std::endl ;
}


Scalar Primal::solve( const SolverOptions& options, DynVec &lambda, SolverStats &stats) const
{
	bogus::Timer timer ;

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

	Scalar res = -1 ;

	if( options.algorithm == SolverOptions::GaussSeidel ) {
		bogus::SOCLaw< 6, Scalar, true > law( m_data.n(), m_data.mu.data() ) ;

		bogus::GaussSeidel< WType > gs( W ) ;
		gs.setTol( options.tolerance );
		gs.setMaxIters( options.maxIterations );
		gs.callback().connect( stats, &SolverStats::ackResidual );

		res = gs.solve( law, m_data.w, lambda ) ;
	} else {
		res = solveCadoux( W, options, lambda, stats ) ;
	}

	stats.residual = res ;
	stats.time = timer.elapsed() ;

	return res ;
}

Scalar Primal::solveCadoux( const WType& W, const SolverOptions& options, DynVec &lambda, SolverStats &stats) const
{
	return -1 ;

}

}
