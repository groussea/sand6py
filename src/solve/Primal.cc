#include "Primal.hh"
#include "PrimalData.hh"

#include "utils/Log.hh"

#include <bogus/Core/Block.impl.hpp>
#include <bogus/Core/BlockSolvers.impl.hpp>
#include <bogus/Extra/SecondOrder.impl.hpp>

#include <bogus/Core/Utils/Timer.hpp>
#include <bogus/Interfaces/Cadoux.hpp>

namespace d6 {

//! Either log directly residual, or re-evaluate it using another complemntarity func
/** (for objective benchmarking purposes) */
template <typename MatrixT >
struct CallbackProxy {

	CallbackProxy( Primal::SolverStats& stats, bogus::Timer &timer,
				   const MatrixT& W, const DynVec& mu, const DynVec &b, const DynVec &x )
		: m_stats( stats ), m_timer( timer ),
		  m_W(W), m_mu(mu), m_b(b), m_x(x)
	{
	}

	void ackResidual( unsigned iter, Scalar err )
	{
		if( m_stats.shouldLogAC ) err = evalAC() ;
		m_stats.log( iter, err, m_timer.elapsed() );
	}

	Scalar evalAC() const {
		const DynVec y = m_W*m_x + m_b ;
		const Index n = m_W.rowsOfBlocks() ;

		Scalar err = 0 ;
#pragma omp parallel for reduction( + : err )
		for( Index i = 0 ; i < n ; ++i ) {
			const VecS lx = m_x.segment<SD>(SD*i) ;
			VecS  ac = lx -   y.segment<SD>(SD*i) ;
			// Normal part
			ac[0] = std::max(0., ac[0])  ;
			//Tangential part
			const Scalar nT = ac.segment<SD-1>(1).norm() ;
			if( nT > lx[0]*m_mu[i] ) {
				ac.segment<SD-1>(1) *= lx[0]*m_mu[i]/nT ;
			}
			//Error
			ac -= lx ;
			const Scalar lerr = ac.squaredNorm() ;
			err += lerr ;
		}
		return err / (1+n) ;
	}

private:
	Primal::SolverStats& m_stats ;
	bogus::Timer& m_timer ;
	bool m_evalAC ;

	const MatrixT & m_W ;
	const DynVec  & m_mu ;
	const DynVec  & m_b ;
	const DynVec  & m_x ;

};

void Primal::SolverStats::log( unsigned iter, Scalar res, Scalar time )
{
	m_log.emplace_back(Primal::SolverStats::Entry{res,time,iter})  ;
	d6::Log::Debug() << "Primal " << iter << " =\t " << res << std::endl ;
	if( time > timeOut ) throw TimeOutException() ;
}



Primal::SolverOptions::SolverOptions()
	: algorithm( GaussSeidel ),
	  maxIterations(250), maxOuterIterations( 15 ),
	  projectedGradientVariant( -1  ),
	  tolerance( 1.e-6 )
{}


Primal::Primal(const PrimalData &data)
	: m_data( data )
{}

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

		bogus::NarySum< JMJtProd > rbSum( m_data.n() * SD, m_data.n() * SD ) ;
		// Add rigid bodies jacobians
		for( unsigned i = 0 ; i < m_data.jacobians.size() ; ++i ) {
			const PrimalData::JacobianType &JM = m_data.jacobians[i] ;
			rbSum +=  JM * ( m_data.inv_inertia_matrices[i] * JM.transpose() ) ;
		}

		typedef bogus::Product< PrimalData::HType, bogus::Transpose< PrimalData::HType > > HHtProd ;
		typedef bogus::Addition< HHtProd, bogus::NarySum< JMJtProd > > Wexpr ;

		const Wexpr W = m_data.H * m_data.H.transpose() + rbSum ;

		bogus::Signal<unsigned, Scalar> callback ;
		CallbackProxy< Wexpr > callbackProxy( stats, timer, W, m_data.mu, m_data.w, lambda ) ;
		callback.connect( callbackProxy, &CallbackProxy<Wexpr>::ackResidual );

		bogus::ProjectedGradient< Wexpr > pg ;

		if( options.projectedGradientVariant < 0 ) {
			pg.setDefaultVariant( bogus::projected_gradient::SPG );
		} else {
			pg.setDefaultVariant( (bogus::projected_gradient::Variant) options.projectedGradientVariant );
		}

		pg.setTol( options.tolerance );
		pg.setMaxIters( options.maxIterations );

		res = bogus::solveCadoux<SD>( W, m_data.w.data(), m_data.mu.data(), pg,
									 lambda.data(), options.maxOuterIterations, &callback ) ;

	} else {

		// Explicit assembly of W

		typedef typename FormMat<SD,SD>::SymType WType ;
		//typedef bogus::SparseBlockMatrix< Eigen::Matrix< Scalar, 6, 6, Eigen::RowMajor >,
	//										bogus::SYMMETRIC > WType ;
		WType W = m_data.H * m_data.H.transpose() ;

		// Add rigid bodies jacobians
		for( unsigned i = 0 ; i < m_data.jacobians.size() ; ++i ) {
			const PrimalData::JacobianType &JM = m_data.jacobians[i] ;

			W +=  JM * m_data.inv_inertia_matrices[i] * JM.transpose() ;
		}

		CallbackProxy< WType > callbackProxy( stats, timer, W, m_data.mu, m_data.w, lambda ) ;

		if( options.algorithm == SolverOptions::GaussSeidel ) {
			// Direct Gauss-Seidel solver

			bogus::SOCLaw< SD, Scalar, true > law( m_data.n(), m_data.mu.data() ) ;

			bogus::GaussSeidel< WType > gs( W ) ;
			gs.setTol( options.tolerance );
			gs.setMaxIters( options.maxIterations );
			gs.callback().connect( callbackProxy, &CallbackProxy<WType>::ackResidual );

			res = gs.solve( law, m_data.w, lambda ) ;

		} else  {

			// Cadoux fixed-point algorithm

			bogus::Signal<unsigned, Scalar> callback ;
			callback.connect( callbackProxy, &CallbackProxy<WType>::ackResidual );

			if( options.algorithm == SolverOptions::Cadoux_GS ) {

				// Gauss-Seidel inner solver

				bogus::GaussSeidel< WType > gs ;
				gs.setTol( options.tolerance );
				gs.setMaxIters( options.maxIterations );

				res = bogus::solveCadoux<SD>( W, m_data.w.data(), m_data.mu.data(), gs,
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

				res = bogus::solveCadoux<SD>( W, m_data.w.data(), m_data.mu.data(), pg,
											  lambda.data(), options.maxOuterIterations, &callback ) ;
			}
		}
	}

	return res ;
}

}
