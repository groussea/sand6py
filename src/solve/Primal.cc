#include "Primal.hh"

#include "utils/Log.hh"

#include <bogus/Core/Block.impl.hpp>
#include <bogus/Core/BlockSolvers.impl.hpp>
#include <bogus/Extra/SecondOrder.impl.hpp>

#include <bogus/Core/Utils/Timer.hpp>

#include <bogus/Interfaces/Cadoux.hpp>

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

	typedef typename FormMat<6,6>::SymType WType ;
	//typedef bogus::SparseBlockMatrix< Eigen::Matrix< Scalar, 6, 6, Eigen::RowMajor >,
//										bogus::SYMMETRIC > WType ;
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
	} else  {

		bogus::Signal<unsigned, Scalar> callback ;
		callback.connect( stats, &Primal::SolverStats::ackResidual );

		if( options.algorithm == SolverOptions::Cadoux_GS ) {

			bogus::GaussSeidel< WType > gs ;
			gs.setTol( options.tolerance );
			gs.setMaxIters( options.maxIterations );

			res = bogus::solveCadoux<6>( W, m_data.w.data(), m_data.mu.data(), gs,
										  lambda.data(), options.maxOuterIterations, &callback ) ;
		} else {

			bogus::ProjectedGradient< WType > pg ;
			pg.setDefaultVariant( bogus::projected_gradient::APGD );
			pg.setTol( options.tolerance );
			pg.setMaxIters( options.maxIterations );

			res = bogus::solveCadoux<6>( W, m_data.w.data(), m_data.mu.data(), pg,
										  lambda.data(), options.maxOuterIterations, &callback ) ;
		}
	}

	stats.residual = res ;
	stats.time = timer.elapsed() ;

	return res ;
}

}
