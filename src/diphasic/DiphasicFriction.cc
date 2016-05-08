#include "DiphasicFriction.hh"


#include <bogus/Core/Block.impl.hpp>
#include <bogus/Core/BlockSolvers.impl.hpp>

#include <bogus/Extra/SecondOrder.impl.hpp>

namespace d6 {


typedef bogus::SparseBlockMatrix<DynMatS> BType ;

static void combine( const DiphasicPrimalData::HType& G,
					 const DiphasicPrimalData::HType& H,
					 const unsigned p,
					 BType& B )
{
	const Index m  = G.colsOfBlocks() ;
	const Index ma = H.colsOfBlocks() ;
	const Index n  = G.rowsOfBlocks() ;

	std::vector< unsigned > cpb ( m+ma+1, WD ) ;
	cpb.back() = p - G.cols() - H.cols()  ;

//	B.setRows( std::vector<unsigned>{ (unsigned)G.rows() } );
//	B.setCols( std::vector<unsigned>{ (unsigned)G.cols(), (unsigned)H.cols(),
//									  (unsigned)(p - G.cols() - H.cols() ) } );
	B.setRows( G.rowsOfBlocks() );
	B.setCols( cpb ) ;

	std::cout << p << " vs " << B.rows() << std::endl ;

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

static void ack( unsigned k, Scalar res ) {
	std::cout << "PG " << k << "\t" << res << std::endl ;
}


Scalar DiphasicFrictionSolver::solve(
		const ESM &M, const DiphasicPrimalData::MInvType& Minv,
		DynVec &x, DynVec &lambda)
{

	typedef DiphasicPrimalData::MInvType MInvType ;

	BType B ;
	combine( m_data.G, m_data.H,  x.rows(), B ) ;


	typedef bogus::Product< BType, bogus::Product< MInvType, bogus::Transpose< BType > > >
			WExpr ;
	WExpr W = B * ( Minv * B.transpose() ) ;

	bogus::SOCLaw< SD, Scalar, true > law( m_data.mu.rows(), m_data.mu.data() ) ;

	bogus::ProjectedGradient< WExpr > pg(W) ;
	pg.callback().connect( &ack );
	pg.useInfinityNorm( true );
	pg.setTol(1.e-12);
	pg.setDefaultVariant( bogus::projected_gradient::SPG );

	const Scalar res = pg.solve( law, m_data.k, lambda) ;
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
