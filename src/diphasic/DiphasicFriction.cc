#include "DiphasicFriction.hh"

#include "utils/Log.hh"

#include <bogus/Core/Block.impl.hpp>
#include <bogus/Core/Block.io.hpp>
#include <bogus/Core/BlockSolvers.impl.hpp>

#include <bogus/Extra/SecondOrder.impl.hpp>
#include <bogus/Core/Utils/Timer.hpp>

#include <bogus/Interfaces/Cadoux.hpp>
#include <bogus/Core/BlockSolvers/ProductGaussSeidel.impl.hpp>

namespace d6 {

DiphasicFrictionSolver::Options::Options()
	: algorithm( PG ),
	  reduced(true), direct(true), useCadoux( true ),
	  maxIterations(250), maxOuterIterations( 15 ),
	  projectedGradientVariant( -1  ),
	  useInfinityNorm( true ), tolerance( 1.e-6 )
{
}

struct DFCallbackProxy {

	DFCallbackProxy( FrictionSolver::Stats& stats, bogus::Timer &timer )
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
static void combineDiag( const DiphasicPrimalData::DType& G,
						 const DiphasicPrimalData::DType& H,
						 const unsigned s,
						 DiphasicPrimalData::DType& B )
{
	const Index m  = G.colsOfBlocks() ;
	const Index n  = H.rowsOfBlocks() ;

	B.setRows( s/WD );
	B.setCols( s/WD ) ;

	for( Index i = 0 ; i < m ; ++i ) {
		for( DiphasicPrimalData::HType::InnerIterator it ( G.majorIndex(), i ) ; it ; ++it  ) {
			B.insertBack( i, it.inner() ) = G.block( it.ptr() ) ;
		}
	}
	for( Index i = 0 ; i < n ; ++i ) {
		for( DiphasicPrimalData::HType::InnerIterator it ( H.majorIndex(), i ) ; it ; ++it  ) {
			B.insertBack( m+i, m+it.inner() ) = H.block( it.ptr() ) ;
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

	DFCallbackProxy callbackProxy( stats, timer ) ;

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
		callback.connect( callbackProxy, &DFCallbackProxy::ackResidual );
		res = bogus::solveCadoux<SD>( W, data.k.data(), data.mu.data(), pg,
									  lambda.data(), options.maxOuterIterations, &callback ) ;
	} else {
		bogus::SOCLaw< SD, Scalar, true > law( data.mu.rows(), data.mu.data() ) ;
		pg.callback().connect( callbackProxy, &DFCallbackProxy::ackResidual );
		res = pg.solve( law, data.k, lambda) ;
	}

	return res ;
}

template <typename WType, typename LSExpr, typename HExpr>
static Scalar solveGS(
		const DiphasicFrictionSolver::Options& options, const DiphasicPrimalData &data,
		const WType& W, const LSExpr& Minv, const HExpr &H,
		DynVec &lambda,
		FrictionSolver::Stats& stats, bogus::Timer& timer	)
{

	DFCallbackProxy callbackProxy( stats, timer ) ;

	bogus::GaussSeidel< WType > gs(W) ;
	gs.useInfinityNorm( options.useInfinityNorm );
	gs.setMaxIters( options.maxIterations );
	gs.setTol( options.tolerance );

	Scalar res = -1 ;

	if( !options.useCadoux ) {
		bogus::SOCLaw< SD, Scalar, true > law( data.mu.rows(), data.mu.data() ) ;
		gs.callback().connect( callbackProxy, &DFCallbackProxy::ackResidual );

		DynVec c = DynVec::Zero( Minv.rows() ) ;
		res = gs.solveWithLinearConstraints( law, Minv, H, 1, data.k, c, lambda, false, 5) ;
	}

	return res ;
}

Scalar DiphasicFrictionSolver::solve(const Options &options,
		DynVec &x, DynVec &lambda, FrictionSolver::Stats &stats ) const
{
	const Scalar pen = 0.e-8 ;

	ESM M ;
	DiphasicPrimalData::MInvType Minv ;


	Scalar res = -1 ;

	if( options.algorithm == Options::ADMM ) {
		res = solveADMM( options, x, lambda, stats ) ;
	} else if( options.reduced ) {
		res = solveRed( options, pen, x, lambda, stats) ;
	} else {

		if(options.direct) {
			m_data.makePenalizedEigenStokesMatrix( M, pen );
			DiphasicPrimalData::factorize( M, Minv ) ;

			res = solveStokes( options, Minv, x, lambda, stats ) ;
		}
	}

	return res ;
}

Scalar DiphasicFrictionSolver::solve(const Options &options,
		const DiphasicPrimalData::MInvType& Minv,
		DynVec &x, DynVec &lambda, FrictionSolver::Stats &stats ) const
{

	if( options.algorithm == Options::ADMM ) {
		return solveADMM( options, x, lambda, stats ) ;
	}
	if( options.direct && !options.reduced ) {
		return solveStokes( options, Minv, x, lambda, stats) ;
	}

	return solve( options, x, lambda, stats )	 ;

}


Scalar DiphasicFrictionSolver::solveStokes(const Options &options,
		const DiphasicPrimalData::MInvType& Minv,
		DynVec &x, DynVec &lambda, FrictionSolver::Stats &stats ) const
{
	assert( options.algorithm == Options::PG && !options.reduced && options.direct ) ;

	bogus::Timer timer ;

	typedef DiphasicPrimalData::MInvType MInvType ;
	typedef DiphasicPrimalData::HType    BType ;

	BType B ;
	combine( m_data.G, m_data.H,  x.rows(), B ) ;


	typedef bogus::Product< BType, bogus::Product< MInvType, bogus::Transpose< BType > > >
			WExpr ;
	WExpr W = B * ( Minv * B.transpose() ) ;


	Scalar res =  -1 ;
	if( options.algorithm == Options::PG ) {
		solvePG( options, m_data, W, lambda, stats, timer ) ;
	}

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

template< typename PInvType, typename KType >
static Scalar solveProdGSRed(const DiphasicFrictionSolver::Options &options,
					   const DiphasicPrimalData& data,
					   const DiphasicPrimalData::DType& R_inv,
					   const PInvType& P_inv, const KType &K,
					   DynVec &lambda, DynVec& delta_p,
					   FrictionSolver::Stats &stats, bogus::Timer &timer )
{
	typedef DiphasicPrimalData::HType HType ;
	typedef DiphasicPrimalData::DType DType ;

	HType H ;
	DType D ;
	combine( data.G, data.H, data.m() + data.r(), H );
	combineDiag( data.M_lumped_inv, R_inv, data.m() + data.r(), D );

	DFCallbackProxy callbackProxy( stats, timer ) ;

	for( Index i = 0 ; i < D.rowsOfBlocks() ; ++i )
		D.block(i).diagonal() = D.block(i).diagonal().array().sqrt() ;
	HType Hsqrt = H * D ;
	bogus::ProductGaussSeidel< HType > gs(Hsqrt) ;

//	bogus::ProductGaussSeidel< HType, DType > gs(H,D) ;

	gs.useInfinityNorm( options.useInfinityNorm );
	gs.setMaxIters( options.maxIterations );
	gs.setTol( options.tolerance );

	Scalar res = -1 ;

	if( !options.useCadoux ) {
		bogus::SOCLaw< SD, Scalar, true > law( data.mu.rows(), data.mu.data() ) ;
		gs.callback().connect( callbackProxy, &DFCallbackProxy::ackResidual );

		DynVec c = DynVec::Zero( P_inv.rows() ) ;
		res = gs.solveWithLinearConstraints( law, P_inv, K, 1, data.k, c, lambda, false, 5) ;
	}

	delta_p = - P_inv * K.transpose() * lambda ;

	return res ;
}

template< typename PInvType, typename KType, typename GAGExpr >
static Scalar solveGSRed(const DiphasicFrictionSolver::Options &options,
					   const DiphasicPrimalData& data,
					   const PInvType& P_inv, const KType &K,
					   const GAGExpr &gag, const GAGExpr& hrh,
					   DynVec &lambda, DynVec& delta_p,
					   FrictionSolver::Stats &stats, bogus::Timer &timer )
{

	typename FormMat<SD,SD>::SymType W = (gag+hrh) ;

	const Scalar res = solveGS( options, data, W, P_inv, K, lambda, stats, timer ) ;

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

	if( options.direct ) {
		typedef DiphasicPrimalData::MInvType PInvType ;

		DiphasicPrimalData::ESM P_csr ;
		bogus::convert( P, P_csr ) ;
		PInvType P_inv ;
		DiphasicPrimalData::factorize( P_csr, P_inv ) ;

		if( options.algorithm == Options::PG ) {
			res = solvePGRed( options, m_data, P_inv, K, gag, hrh, lambda, delta_p, stats, timer ) ;
		} else if( options.algorithm == Options::GS ) {
//			res = solveGSRed( options, m_data, P_inv, K, gag, hrh, lambda, delta_p, stats, timer ) ;
			res = solveProdGSRed( options, m_data, R_inv, P_inv, K, lambda, delta_p, stats, timer ) ;
		}
	} else {
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

		if( options.algorithm == Options::PG ) {
			res = solvePGRed( options, m_data, P_inv, K, gag, hrh, lambda, delta_p, stats, timer ) ;
		} else if( options.algorithm == Options::GS ) {
			res = solveGSRed( options, m_data, P_inv, K, gag, hrh, lambda, delta_p, stats, timer ) ;
		}
	}


	x.segment( m_data.m() + m_data.r(), m_data.p() ) += delta_p ;
	x.head( m_data.m() ) += m_data.M_lumped_inv *
			( m_data.B.transpose() * delta_p + m_data.G.transpose() * lambda) ;
	x.segment( m_data.m(), m_data.r() ) += R_inv *
			( m_data.C.transpose() * delta_p + m_data.H.transpose() * lambda) ;


	return res ;
}


Scalar DiphasicFrictionSolver::solveADMM(const Options &options,
		DynVec &x, DynVec &lambda, FrictionSolver::Stats &stats ) const
{
	bogus::Timer timer ;
	DFCallbackProxy callbackProxy( stats, timer ) ;

	Scalar res = -1 ;

	//// Combinations ////
	bogus::SparseBlockMatrix< DiphasicPrimalData::AType > A ;
	A.setRows( std::vector<unsigned>{ (unsigned)m_data.m(), (unsigned)m_data.r() });
	A.setCols( std::vector<unsigned>{ (unsigned)m_data.m(), (unsigned)m_data.r() });
	A.insertBack(0,0) = m_data.A ;
	A.insertBack(1,1) = m_data.R ;
	A.finalize();

	bogus::SparseBlockMatrix< DiphasicPrimalData::CType > B ;
	B.setRows( std::vector<unsigned>{ (unsigned)m_data.B.rows() });
	B.setCols( std::vector<unsigned>{ (unsigned)m_data.m(), (unsigned)m_data.r() });
	B.insertBack(0,0) = m_data.B ;
	B.insertBack(0,1) = m_data.C ;
	B.finalize();

	bogus::SparseBlockMatrix< DiphasicPrimalData::HType > H ;
	H.setRows( std::vector<unsigned>{ (unsigned)m_data.G.rows() });
	H.setCols( std::vector<unsigned>{ (unsigned)m_data.m(), (unsigned)m_data.r() });
	H.insertBack(0,0) = m_data.G ;
	H.insertBack(0,1) = m_data.H ;
	H.finalize();


	DiphasicPrimalData::DType P ;
	P.setRows( m_data.A.rowsOfBlocks() + m_data.R.rowsOfBlocks() ) ;
	P.setIdentity() ;

	for( Index i = 0 ; i < P.rowsOfBlocks() ; ++i ) {
		Index j = i - m_data.A.rowsOfBlocks();
		if( j < 0 ) {
			P.diagonal(i).diagonal() = m_data.M_lumped_inv.diagonal(i).diagonal() ;
		} else {
			P.block(j).diagonal() = 1./(1.e-8 + m_data.R.block(j).diagonal().array()) ;
		}
	}

	/////////////////////

//	Scalar max_R = 0 ;
//	for( Index i = 0 ;  i < m_data.R.nBlocks() ; ++i ) {
//		max_R = std::max( max_R, m_data.R.block(i).lpNorm<Eigen::Infinity>() ) ;
//	}
//	Scalar max_A = 0 ;
//	for( Index i = 0 ;  i < m_data.A.nBlocks() ; ++i ) {
//		max_A = std::max( max_A, m_data.A.block(i).lpNorm<Eigen::Infinity>() ) ;
//	}
//	std::cout << "R " << max_R << "  A " << max_A << std::endl ;

	/////////////////////

	const DynVec f = DynVec::Zero( A.rows() ) ;
	const DynVec b = DynVec::Zero( B.rows() ) ;

	DynVec v = DynVec::Zero( A.rows() ) ; //FIXME warm start
	DynVec p = DynVec::Zero( B.rows() ) ;


	bogus::DualAMA< DiphasicPrimalData::HType > dama( m_data.G ) ;
	dama.useInfinityNorm( options.useInfinityNorm );
	dama.setMaxIters( options.maxIterations );
	dama.setTol( options.tolerance );

	bogus::SOCLaw< SD, Scalar, true > law( m_data.mu.rows(), m_data.mu.data() ) ;
	dama.callback().connect( callbackProxy, &DFCallbackProxy::ackResidual );


//	dama.setDefaultVariant( bogus::admm::Accelerated );
	dama.setLineSearchIterations( 0 );
//	dama.setFpStepSize( 0.3 );
//	dama.setProjStepSize( 2 );
	dama.setFpStepSize( 5.e-2 );
	dama.setProjStepSize( 30. );

	typedef bogus::MatrixPreconditioner< DiphasicPrimalData::DType >
			::Type < bogus::BlockObjectBase<bogus::SparseBlockMatrix< DiphasicPrimalData::AType >> >
			PrecondType ;
	PrecondType precond ;
	precond.setPreconditionerMatrix( P );
//	typedef bogus::TrivialPreconditioner<
//			bogus::BlockObjectBase<bogus::SparseBlockMatrix< DiphasicPrimalData::AType >> >
//			PrecondType ;
//	PrecondType precond ;

	DynVec test_0 ( A.rows() ), test_1( A.rows() ) ;
	test_0.setOnes() ;
	precond.apply< false >( test_0, test_1 ) ;
	std::cout << "C " << test_1.minCoeff() << "  " << test_1.maxCoeff() << std::endl ;

	res = dama.solveWithLinearConstraints< bogus::admm::Standard>
			( law, A, B, H, precond, f, b, m_data.k, v, p, lambda, .003 ) ;

	x.head( A.rows() ) += v ;
	x.segment( A.rows(), B.rows() ) += p ;

	return res ;
}


} //d6
