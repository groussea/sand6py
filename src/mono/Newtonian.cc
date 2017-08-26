#include "Newtonian.hh"
#include "PhaseStepData.hh"


#include "utils/Log.hh"
#include "utils/Config.hh"
#include "utils/Stats.hh"

#include "simu/FormBuilder.hh"

#include <bogus/Core/Block.impl.hpp>
#include <bogus/Core/Block.io.hpp>
#include <bogus/Core/BlockSolvers/Krylov.impl.hpp>
#include <bogus/Core/Utils/Timer.hpp>

// phi Div u = 0    ->  int grad( phi p ). v

// Ghost fluid, alpha = +inf
//     Div u = 0    ->  int ( phi grad( p ) + (1-phi) grad(p) ). v
//                   =  int (grad p)


namespace d6 {


namespace
{

typedef FormMat<WD,WD>::Type AType ;
typedef FormMat< 1,WD>::Type BType ;

typedef Eigen::SparseMatrix< Scalar > ESM ;
typedef Eigen::SimplicialLDLT<ESM> ESMFact ;
typedef bogus::SparseBlockMatrix< bogus::EigenSparseFactorization< ESM, ESMFact > > MInvType ;

void makePenalizedEigenStokesMatrix(
        const AType& iA,
        const BType& iB,
        const Scalar pen,
        ESM& M
        )
{
	// Assemble Eigen mat for
	// (  A  -B^T  )
	// ( -B  -pen Id  )

	const Index m = iA.rows() ;
	const Index p = iB.rows() ;
	const Index s = m+p ;

	typedef Eigen::SparseMatrix< Scalar > SM ;
	SM A, B ;
	bogus::convert( iA, A ) ;
	bogus::convert( iB, B ) ;

	A.prune(1.) ;
	B.prune(1.) ;

	M.resize( s, s ) ;

	typedef Eigen::Triplet<Scalar> Tpl ;
	std::vector< Tpl > tpl ;
	tpl.reserve( A.nonZeros() + p + 2*B.nonZeros()  );

	for( Index i = 0 ; i < m ; ++i ) {
		for( SM::InnerIterator it (A, i) ; it ; ++it )
		{
			tpl.push_back( Tpl( it.row(), i, it.value() ) );
		}
		for( SM::InnerIterator it (B, i) ; it ; ++it )
		{
			tpl.push_back( Tpl( m + it.row(), i, -it.value() ) );
			tpl.push_back( Tpl( i, m + it.row(), -it.value() ) );
		}
	}
	for( Index i = 0 ; i < p ; ++i ) {
//		for( SM::InnerIterator it (Q, i) ; it ; ++it )
//		{
//			tpl.push_back( Tpl( m+r+it.row(), m+r+i, -pen * it.value() ) );
//		}
		tpl.push_back( Tpl( m + i, m + i, - pen ) );
	}

	M.setFromTriplets( tpl.begin(), tpl.end() ) ;
}

void  factorize( ESM &mat, MInvType &sbm )
{

	std::vector< unsigned > rpb {(unsigned)  mat.rows()} ;
	std::vector< unsigned > cpb {(unsigned)  mat.cols()} ;

	sbm.setRows( rpb );
	sbm.setCols( cpb );

	sbm.insertBack(0,0).compute( mat ) ;
	sbm.finalize();

}


}


void NewtonianSolver::solveIncompressibility(
        const Config &c, const Scalar, const PhaseStepData& stepData,
        std::vector<RigidBodyData> &,
        DynVec &u, Phase& , Stats &simuStats )
{
	// solve A delta_u - B' p = 0
	//       - delta_u        = B u

	typedef FormMat<1,WD>::Type BType ;
	BType  B = stepData.proj.pressure * ( stepData.forms.B * stepData.proj.vel ) ;
	DynVec w = stepData.proj.pressure * ( stepData.forms.B * u )  ;

	const Index m  = B.cols() ;
	const Index p  = B.rows() ;

	DynVec x ( m+p );
	x.setZero() ;
	x.segment(m,p).setZero() ; // = fluid.pressure.flatten() ;

	DynVec l (m+p ) ;
	l.setZero() ;
	l.segment(m,p) = w ;

	bogus::Timer timer ;

	bool iterative_stokes = false ;
	if(iterative_stokes)
	{
		typedef FormMat<WD,WD>::Type MInvType ;
		const MInvType &M_fac = stepData.forms.M_lumped_inv  ;
		// solve
		//   Mu - B'p = lu
		//  -Bu - c p = lp
		//
		//  cp + B ( M^-1 lu + M^-1 B' p ) + lp = 0
		//  W p - b = 0

		typedef bogus::Product<
		        bogus::Product < BType, MInvType >,
		        bogus::Transpose< BType > > WType ;

		WType W = B * M_fac * B.transpose() ;
		Eigen::VectorXd b = - ( B * l.head(m) + l.segment(m,p) ) ;

		typedef Eigen::Matrix<Scalar,1,1> ScalarAsMat ;
		typedef bogus::SparseBlockMatrix< ScalarAsMat, bogus::SYMMETRIC > PenType ;
		PenType pen ; pen.setRows( W.rows() ) ; pen.setIdentity() ;

		typedef bogus::Addition< WType, bogus::Scaling<PenType> > WPenType ;
		WPenType WPen = W + c.compressibility * pen ;

		bogus::Krylov< WPenType > krylov( WPen ) ;
		krylov.setMaxIters(100);
		krylov.setTol(1.e-8);
		auto xp = x.segment( m, p ) ;
		double res = krylov.asCG().solve( b, xp ) ;

		x.head( m ) = M_fac * ( l.head(m) + B.transpose()*xp ) ;
		Log::Debug() << "Krylov Stokes solver time: " << timer.elapsed() <<" \t res: " << res << std::endl ;

	} else {
		ESM M ;
		makePenalizedEigenStokesMatrix(stepData.forms.A, B, c.compressibility, M);
		MInvType M_fac ;
		factorize(M, M_fac);
		Log::Debug() << "Stokes factorization time: " << timer.elapsed() << std::endl ;

		if( M_fac.block(0).factorization().info() != 0 ) {
			Log::Error() << "Stokes fac failed! "  << std::endl ;
			std::abort() ;
		}

		x = M_fac * l ;
		Log::Debug() << "Direct Stokes solver time: " << timer.elapsed() << std::endl ;
	}

	simuStats.lcpSolveTime = timer.elapsed() ;

	// Update velocity
	u += x.head(m)  ;
}


}
