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
static Scalar solvePG(
		const DiphasicFrictionSolver::Options& options, const DiphasicPrimalData &data,
		const WExpr& W, DynVec &lambda,
		FrictionSolver::Stats& stats, bogus::Timer& timer	)
{

	CallbackProxy callbackProxy( stats, timer ) ;

	bogus::ProjectedGradient< WExpr > pg(W) ;
	pg.useInfinityNorm( options.useInfinityNorm );
	pg.setMaxIters( options.maxIterations );
	pg.setTol( options.tolerance );

	if( options.projectedGradientVariant < 0 ) {
		pg.setDefaultVariant( bogus::projected_gradient::SPG );
	} else {
		pg.setDefaultVariant( (bogus::projected_gradient::Variant) options.projectedGradientVariant );
	}

	Scalar res = -1 ;

	if( options.useCadoux ) {
		bogus::Signal<unsigned, Scalar> callback ;
		callback.connect( callbackProxy, &CallbackProxy::ackResidual );
		res = bogus::solveCadoux<SD>( W, data.k.data(), data.mu.data(), pg,
									  lambda.data(), options.maxOuterIterations, &callback ) ;
	} else {
		bogus::SOCLaw< SD, Scalar, true > law( data.mu.rows(), data.mu.data() ) ;
		pg.callback().connect( callbackProxy, &CallbackProxy::ackResidual );
		res = pg.solve( law, data.k, lambda) ;
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

		res = solveStokes( options, Minv, x, lambda, stats ) ;
	}
	if( options.algorithm == Options::PG_Fac_Red || options.algorithm == Options::PG_CG_Red ) {
		res = solveRed( options, pen, x, lambda, stats) ;
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


	Scalar res = solvePG( options, m_data, W, lambda, stats, timer ) ;

	x += Minv * B.transpose() * lambda ;

	return res ;

}

template< typename PInvType, typename KType, typename GAGExpr >
static Scalar solvePGRed(const DiphasicFrictionSolver::Options &options,
					   const DiphasicPrimalData& data,
					   const PInvType& P_inv, const KType &K,
					   const GAGExpr &gag, const GAGExpr& hrh,
					   DynVec &lambda, DynVec& delta_p,
					   FrictionSolver::Stats &stats, bogus::Timer &timer )
{


	typedef bogus::Product< KType, bogus::Product< PInvType, bogus::Transpose< KType > > >
			BPBExpr ;

	BPBExpr kpk = K * ( P_inv * K.transpose() ) ;

	typedef bogus::Addition< bogus::Addition< GAGExpr, GAGExpr >, BPBExpr > WExpr ;
//	typedef bogus::Addition< GAGExpr, GAGExpr > WExpr ;
	WExpr W = (gag+hrh) - kpk ;

	const Scalar res = solvePG( options, data, W, lambda, stats, timer ) ;

	delta_p = - P_inv * K.transpose() * lambda ;

	return res ;
}


Scalar DiphasicFrictionSolver::solveRed(const Options &options, const Scalar pen,
		DynVec &x, DynVec &lambda, FrictionSolver::Stats &stats ) const
{

	bogus::Timer timer ;
	Scalar res = -1 ;

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

	BType K  = m_data.G * ( m_data.M_lumped_inv * m_data.B.transpose() ) +
			   m_data.H * (        R_inv        * m_data.C.transpose() ) ;

	// TODO use SymType and SelfAdjointView
	typedef FormMat<1,1>::Type PType ;
	PType P =  m_data.B * ( m_data.M_lumped_inv * m_data.B.transpose() ) +
			   m_data.C * (        R_inv        * m_data.C.transpose() ) ;
	P  += pen * P.Identity() ;

	DynVec delta_p ;

	if( options.algorithm == Options::PG_Fac_Red ) {
		typedef DiphasicPrimalData::MInvType PInvType ;

		DiphasicPrimalData::ESM P_csr ;
		bogus::convert( P, P_csr ) ;
		PInvType P_inv ;
		DiphasicPrimalData::factorize( P_csr, P_inv ) ;

		res = solvePGRed( options, m_data, P_inv, K, gag, hrh, lambda, delta_p, stats, timer ) ;
	} else if( options.algorithm == Options::PG_CG_Red ) {
		typedef bogus::Krylov< PType >::CGType CGType ;
		typedef bogus::SparseBlockMatrix< CGType > PInvType ;

		bogus::Krylov< PType > cg( P ) ;
		cg.setTol( 1.e-8 );
		cg.setMaxIters( 100 );

		PInvType P_inv ;
		P_inv.setRows(std::vector<unsigned>{(unsigned)P.rows()});
		P_inv.setCols(std::vector<unsigned>{(unsigned)P.rows()});
		P_inv.insertBack(0,0) = cg.asCG().enableResCaching() ;
		P_inv.finalize();

		res = solvePGRed( options, m_data, P_inv, K, gag, hrh, lambda, delta_p, stats, timer ) ;
	}


	x.segment( m_data.m() + m_data.r(), m_data.p() ) += delta_p ;
	x.head( m_data.m() ) += m_data.M_lumped_inv *
			( m_data.B.transpose() * delta_p + m_data.G.transpose() * lambda) ;
	x.segment( m_data.m(), m_data.r() ) += R_inv *
			( m_data.C.transpose() * delta_p + m_data.H.transpose() * lambda) ;


	return res ;
}


} //d6
