#include "Primal.hh"
#include "PrimalData.hh"

#include "utils/Log.hh"

#include <bogus/Core/Block.impl.hpp>
#include <bogus/Core/BlockSolvers.impl.hpp>
#include <bogus/Extra/SecondOrder.impl.hpp>

#include <bogus/Core/Utils/Timer.hpp>
#include <bogus/Interfaces/Cadoux.hpp>

namespace d6 {


Primal::SolverOptions::SolverOptions()
	: algorithm( GaussSeidel ),
	  maxIterations(250), maxOuterIterations( 15 ),
	  projectedGradientVariant( -1  ),
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
//	m_data.dump("last.d6") ;

	bogus::Timer timer ;
	Scalar res = -1 ;

	if( options.algorithm == SolverOptions::Cadoux_PG_NoAssembly ) {

		// Keep W as an expression
		typedef bogus::Product< PrimalData::JacobianType,
				bogus::Product< PrimalData::InvInertiaType, bogus::Transpose< PrimalData::JacobianType > > >
				JMJtProd ;

		bogus::NarySum< JMJtProd > rbSum( m_data.n() * 6, m_data.n() * 6 ) ;
		// Add rigid bodies jacobians
		for( unsigned i = 0 ; i < m_data.jacobians.size() ; ++i ) {
			const PrimalData::JacobianType &JM = m_data.jacobians[i] ;
			rbSum +=  JM * ( m_data.inv_inertia_matrices[i] * JM.transpose() ) ;
		}

		typedef bogus::Product< PrimalData::HType, bogus::Transpose< PrimalData::HType > > HHtProd ;
		typedef bogus::Addition< HHtProd, bogus::NarySum< JMJtProd > > Wexpr ;

		const Wexpr W = m_data.H * m_data.H.transpose() + rbSum ;

		bogus::Signal<unsigned, Scalar> callback ;
		callback.connect( stats, &Primal::SolverStats::ackResidual );
		bogus::ProjectedGradient< Wexpr > pg ;

		if( options.projectedGradientVariant < 0 ) {
			pg.setDefaultVariant( bogus::projected_gradient::SPG );
		} else {
			pg.setDefaultVariant( (bogus::projected_gradient::Variant) options.projectedGradientVariant );
		}

		pg.setTol( options.tolerance );
		pg.setMaxIters( options.maxIterations );

		res = bogus::solveCadoux<6>( W, m_data.w.data(), m_data.mu.data(), pg,
									 lambda.data(), options.maxOuterIterations, &callback ) ;

	} else {

		// Explicit assembly of W

		typedef typename FormMat<6,6>::SymType WType ;
		//typedef bogus::SparseBlockMatrix< Eigen::Matrix< Scalar, 6, 6, Eigen::RowMajor >,
	//										bogus::SYMMETRIC > WType ;
		WType W = m_data.H * m_data.H.transpose() ;

		// Add rigid bodies jacobians
		for( unsigned i = 0 ; i < m_data.jacobians.size() ; ++i ) {
			const PrimalData::JacobianType &JM = m_data.jacobians[i] ;

			W +=  JM * m_data.inv_inertia_matrices[i] * JM.transpose() ;
		}

		if( options.algorithm == SolverOptions::GaussSeidel ) {
			// Direct Gauss-Seidel solver

			bogus::SOCLaw< 6, Scalar, true > law( m_data.n(), m_data.mu.data() ) ;

			bogus::GaussSeidel< WType > gs( W ) ;
			gs.setTol( options.tolerance );
			gs.setMaxIters( options.maxIterations );
			gs.callback().connect( stats, &SolverStats::ackResidual );

			res = gs.solve( law, m_data.w, lambda ) ;

		} else  {

			// Cadoux fixed-point algorithm

			bogus::Signal<unsigned, Scalar> callback ;
			callback.connect( stats, &Primal::SolverStats::ackResidual );

			if( options.algorithm == SolverOptions::Cadoux_GS ) {

				// Gauss-Seidel inner solver

				bogus::GaussSeidel< WType > gs ;
				gs.setTol( options.tolerance );
				gs.setMaxIters( options.maxIterations );

				res = bogus::solveCadoux<6>( W, m_data.w.data(), m_data.mu.data(), gs,
											  lambda.data(), options.maxOuterIterations, &callback ) ;
			} else {

				// Projected Gradient inner solver

				bogus::ProjectedGradient< WType > pg ;

				if( options.projectedGradientVariant < 0 ) {
					pg.setDefaultVariant( bogus::projected_gradient::SPG );
				} else {
					pg.setDefaultVariant( (bogus::projected_gradient::Variant) options.projectedGradientVariant );
				}

				pg.setTol( options.tolerance );
				pg.setMaxIters( options.maxIterations );

				res = bogus::solveCadoux<6>( W, m_data.w.data(), m_data.mu.data(), pg,
											  lambda.data(), options.maxOuterIterations, &callback ) ;
			}
		}
	}

	stats.residual = res ;
	stats.time = timer.elapsed() ;

	return res ;
}

}
