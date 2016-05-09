#include "DiphasicFriction.hh"

#include "utils/Log.hh"

#include <bogus/Core/Block.impl.hpp>
#include <bogus/Core/BlockSolvers.impl.hpp>

#include <bogus/Extra/SecondOrder.impl.hpp>
#include <bogus/Core/Utils/Timer.hpp>

#include <bogus/Interfaces/Cadoux.hpp>

namespace d6 {

DiphasicFrictionSolver::Options::Options()
	: algorithm( PG_Fac_Stokes ), useCadoux( true ) ,
	  maxIterations(250), maxOuterIterations( 15 ),
	  projectedGradientVariant( -1  ),
	  useInfinityNorm( true ), tolerance( 1.e-6 )
{
}

struct CallbackProxy {

	CallbackProxy( FrictionSolver::Stats& stats, bogus::Timer &timer )
		: m_stats( stats ), m_timer( timer )
	{
	}

	void ackResidual( unsigned iter, Scalar err )
	{
		m_stats.log( iter, err, m_timer.elapsed() );
	}

private:
	FrictionSolver::Stats& m_stats ;
	bogus::Timer& m_timer ;

} ;

static void combine( const DiphasicPrimalData::HType& G,
					 const DiphasicPrimalData::HType& H,
					 const unsigned s,
					 DiphasicPrimalData::HType& B )
{
	const Index m  = G.colsOfBlocks() ;
	const Index n  = G.rowsOfBlocks() ;

	B.setRows( n );
	B.setCols( s/WD ) ;

	for( Index i = 0 ; i < n ; ++i ) {
		for( DiphasicPrimalData::HType::InnerIterator it ( G.majorIndex(), i ) ; it ; ++it  ) {
			B.insertBack( i, it.inner() ) = G.block( it.ptr() ) ;
		}
		for( DiphasicPrimalData::HType::InnerIterator it ( H.majorIndex(), i ) ; it ; ++it  ) {
			B.insertBack( i, m+it.inner() ) = H.block( it.ptr() ) ;
		}
	}
	B.finalize();
}

Scalar DiphasicFrictionSolver::solve(const Options &options,
		DynVec &x, DynVec &lambda, FrictionSolver::Stats &stats )
{
	ESM M ;
	DiphasicPrimalData::MInvType Minv ;

	if( options.algorithm == Options::PG_Fac_Stokes) {
		m_data.makePenalizedEigenStokesMatrix( M, 1.e-8 );
		DiphasicPrimalData::factorize( M, Minv ) ;
	}

	return solve( options, M, Minv, x, lambda, stats ) ;

}

Scalar DiphasicFrictionSolver::solve(const Options &options,
		const ESM &M, const DiphasicPrimalData::MInvType& Minv,
		DynVec &x, DynVec &lambda, FrictionSolver::Stats &stats )
{
	if( options.algorithm != Options::PG_Fac_Stokes ) {
		Log::Error() << "Not impl" << std::endl ;
		std::exit(1) ;
	}

	bogus::Timer timer ;

	typedef DiphasicPrimalData::MInvType MInvType ;
	typedef DiphasicPrimalData::HType    BType ;

	BType B ;
	combine( m_data.G, m_data.H,  x.rows(), B ) ;


	typedef bogus::Product< BType, bogus::Product< MInvType, bogus::Transpose< BType > > >
			WExpr ;
	WExpr W = B * ( Minv * B.transpose() ) ;

	bogus::SOCLaw< SD, Scalar, true > law( m_data.mu.rows(), m_data.mu.data() ) ;

	CallbackProxy callbackProxy( stats, timer ) ;

	bogus::ProjectedGradient< WExpr > pg(W) ;
	pg.useInfinityNorm( options.useInfinityNorm );
	pg.setMaxIters( options.maxIterations );
	pg.setTol( options.tolerance );

	if( options.projectedGradientVariant < 0 ) {
		pg.setDefaultVariant( bogus::projected_gradient::Conjugated );
	} else {
		pg.setDefaultVariant( (bogus::projected_gradient::Variant) options.projectedGradientVariant );
	}

	Scalar res = -1 ;

	if( options.useCadoux ) {
		bogus::Signal<unsigned, Scalar> callback ;
		callback.connect( callbackProxy, &CallbackProxy::ackResidual );
		res = bogus::solveCadoux<SD>( W, m_data.k.data(), m_data.mu.data(), pg,
									  lambda.data(), options.maxOuterIterations, &callback ) ;
	} else {
		pg.callback().connect( callbackProxy, &CallbackProxy::ackResidual );
		res = pg.solve( law, m_data.k, lambda) ;
	}

	std::cout << "Res " << res << std::endl ;
	std::cout << "LLLL " << lambda.lpNorm<Eigen::Infinity>() << std::endl ;


	DynVec GHx  = m_data.G*x.head( m_data.G.cols() ) + m_data.H*x.segment( m_data.G.cols(),m_data.H.cols()) ;
	DynVec   Bx  = B*x ;
	std::cout << "TESTB " << (GHx - Bx).squaredNorm() << std::endl ;
	DynVec Bl  = B.transpose()*lambda ;
	DynVec MiBl  = Minv*Bl ;
	DynVec MMiBl = M*MiBl ;
	std::cout << "TESTM " << (MMiBl - Bl).squaredNorm() << std::endl ;

	x += Minv * B.transpose() * lambda ;

	return res ;

}


} //d6
