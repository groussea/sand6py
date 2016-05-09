#include "DiphasicFriction.hh"

#include "utils/Log.hh"

#include <bogus/Core/Block.impl.hpp>
#include <bogus/Core/Block.io.hpp>
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

template <typename WExpr>
Scalar DiphasicFrictionSolver::solvePG( const Options& options, const WExpr& W, DynVec &lambda,
				FrictionSolver::Stats& stats, bogus::Timer& timer	) const
{
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

	return res ;
}

Scalar DiphasicFrictionSolver::solve(const Options &options,
		DynVec &x, DynVec &lambda, FrictionSolver::Stats &stats ) const
{
	const Scalar pen = 1.e-8 ;

	ESM M ;
	DiphasicPrimalData::MInvType Minv ;


	Scalar res = -1 ;

	if( options.algorithm == Options::PG_Fac_Stokes) {
		m_data.makePenalizedEigenStokesMatrix( M, pen );
		DiphasicPrimalData::factorize( M, Minv ) ;

		solveStokes( options, Minv, x, lambda, stats ) ;
	}
	if( options.algorithm == Options::PG_Fac_Red) {
		solveRed( options, pen, x, lambda, stats) ;
	}



	return res ;
}

Scalar DiphasicFrictionSolver::solve(const Options &options,
		const ESM &, const DiphasicPrimalData::MInvType& Minv,
		DynVec &x, DynVec &lambda, FrictionSolver::Stats &stats ) const
{

	if( options.algorithm == Options::PG_Fac_Stokes ) {
		return solveStokes( options, Minv, x, lambda, stats) ;
	}

	return solve( options, x, lambda, stats )	 ;

}


Scalar DiphasicFrictionSolver::solveStokes(const Options &options,
		const DiphasicPrimalData::MInvType& Minv,
		DynVec &x, DynVec &lambda, FrictionSolver::Stats &stats ) const
{
	assert( options.algorithm == Options::PG_Fac_Stokes ) ;

	bogus::Timer timer ;

	typedef DiphasicPrimalData::MInvType MInvType ;
	typedef DiphasicPrimalData::HType    BType ;

	BType B ;
	combine( m_data.G, m_data.H,  x.rows(), B ) ;


	typedef bogus::Product< BType, bogus::Product< MInvType, bogus::Transpose< BType > > >
			WExpr ;
	WExpr W = B * ( Minv * B.transpose() ) ;


	Scalar res = solvePG( options, W, lambda, stats, timer ) ;

	x += Minv * B.transpose() * lambda ;

	return res ;

}

Scalar DiphasicFrictionSolver::solveRed(const Options &options, const Scalar pen,
		DynVec &x, DynVec &lambda, FrictionSolver::Stats &stats ) const
{

	bogus::Timer timer ;

	// R_inv
	DiphasicPrimalData::DType R_inv = m_data.R.Identity() ;
	for( unsigned i = 0 ; i < R_inv.nBlocks() ; ++i )	{
		R_inv.block(i).diagonal().array() = 1./( m_data.R.diagonal(i).diagonal().array() + pen ) ;
	}

	typedef DiphasicPrimalData::HType    HType ;
	typedef bogus::Product< HType, bogus::Product< DiphasicPrimalData::DType, bogus::Transpose< HType > > >
			GAGExpr ;
	GAGExpr gag = m_data.G * ( m_data.M_lumped_inv * m_data.G.transpose() ) ;
	GAGExpr hrh = m_data.H * (        R_inv        * m_data.H.transpose() ) ;


	typedef FormMat< SD, 1>::Type        BType ;
	typedef DiphasicPrimalData::MInvType MInvType ;

	BType K  = m_data.G * ( m_data.M_lumped_inv * m_data.B.transpose() ) ;
		  K += m_data.H * (        R_inv        * m_data.C.transpose() ) ;

	DiphasicPrimalData::ESM P ;
	{
		// TODO use SymType and SelfAdjointView
		FormMat<1,1>::Type P_bsr =  m_data.B * ( m_data.M_lumped_inv * m_data.B.transpose() ) +
									m_data.C * (        R_inv        * m_data.C.transpose() ) ;
		P_bsr += pen * P_bsr.Identity() ;
		bogus::convert( P_bsr, P ) ;
	}


	MInvType P_inv ;
	DiphasicPrimalData::factorize( P, P_inv ) ;

	typedef bogus::Product< BType, bogus::Product< MInvType, bogus::Transpose< BType > > >
			BPBExpr ;

	BPBExpr kpk = K * ( P_inv * K.transpose() ) ;

	typedef bogus::Addition< bogus::Addition< GAGExpr, GAGExpr >, BPBExpr > WExpr ;
//	typedef bogus::Addition< GAGExpr, GAGExpr > WExpr ;
	WExpr W = (gag+hrh) - kpk ;


	Scalar res = solvePG( options, W, lambda, stats, timer ) ;


	const DynVec delta_p = - P_inv * K.transpose() * lambda ;

	x.segment( m_data.m() + m_data.r(), m_data.p() ) += delta_p ;
	x.head( m_data.m() ) += m_data.M_lumped_inv *
			( m_data.B.transpose() * delta_p + m_data.G.transpose() * lambda) ;
	x.segment( m_data.m(), m_data.r() ) += R_inv *
			( m_data.C.transpose() * delta_p + m_data.H.transpose() * lambda) ;


	return res ;
}


} //d6
