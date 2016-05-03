#include "DiphasicSolver.hh"

#include "DiphasicStepData.hh"

#include "mono/Phase.hh"

#include "simu/LinearSolver.hh"

#include "utils/Config.hh"
#include "utils/Log.hh"

#include <bogus/Core/Utils/Timer.hpp>
#include <bogus/Core/Block.impl.hpp>
#include <bogus/Core/Block.io.hpp>


#include <Eigen/Sparse>

namespace d6 {


static void makePenalizedEigenStokesMatrix(
		const DiphasicStepData& stepData, Eigen::SparseMatrix<Scalar> & M,
		const Scalar pen = 1.e-8
		)
{
	// Assemble Eigen mat for
	// (  A    0     -B^T  )
	// (  0    R     -C^T  )
	// ( -B   -C  -pen Id  )

	const Index m  = stepData.forms.A.rows() ;
	const Index ma = stepData.forms.R.rows() ;
	const Index n  = stepData.forms.B.rows() ;
	const Index r = m + n + ma;

	typedef Eigen::SparseMatrix< Scalar > SM ;
	SM A, B, C, R ;
	bogus::convert( stepData.forms.A, A ) ;
	bogus::convert( stepData.forms.R, R ) ;

	//FIXME bogus single-line assignment
	{
		typename FormMat<1,WD>::Type Bproj ;
		Bproj = stepData.forms.B * stepData.fullGridProj.vel ;
		bogus::convert( Bproj, B ) ;

		typename FormMat<1,WD>::Type Cproj ;
		Cproj = stepData.forms.C * stepData.activeProj.vel ;
		bogus::convert( Cproj, C ) ;
	}

	A.prune(1.) ;
	B.prune(1.) ;
	C.prune(1.) ;
	R.prune(1.) ;

	M.resize( r, r ) ;

	typedef Eigen::Triplet<Scalar> Tpl ;
	std::vector< Tpl > tpl ;
	tpl.reserve( A.nonZeros() + R.nonZeros() + n + 2*C.nonZeros() + 2*B.nonZeros()  );

	for( Index i = 0 ; i < m ; ++i ) {
		for( SM::InnerIterator it (A, i) ; it ; ++it )
		{
			tpl.push_back( Tpl( it.row(), i, it.value() ) );
		}
		for( SM::InnerIterator it (B, i) ; it ; ++it )
		{
			tpl.push_back( Tpl( m+ma + it.row(), i, -it.value() ) );
			tpl.push_back( Tpl( i, m+ma + it.row(), -it.value() ) );
		}
	}
	for( Index i = 0 ; i < ma ; ++i ) {
		for( SM::InnerIterator it (R, i) ; it ; ++it )
		{
			tpl.push_back( Tpl( m+it.row(), m+i, it.value() ) );
		}
		for( SM::InnerIterator it (C, i) ; it ; ++it )
		{
			tpl.push_back( Tpl( m+ma + it.row(), m+i, -it.value() ) );
			tpl.push_back( Tpl( m+i, m+ma + it.row(), -it.value() ) );
		}
	}
	for( Index i = 0 ; i < n ; ++i ) {
		tpl.push_back( Tpl( m+ma + i, m+ma + i, - pen ) );
	}

	M.setFromTriplets( tpl.begin(), tpl.end() ) ;

}



DiphasicSolver::DiphasicSolver(const DynParticles &particles)
	: m_particles(particles)
{

}

void DiphasicSolver::step(const Config &config, const Scalar dt, Phase &phase ) const
{

	DiphasicStepData stepData ;
	stepData.compute( m_particles, config, dt, phase );

	solve( config, dt, stepData, phase) ;
}


void DiphasicSolver::solve(
	const Config& config, const Scalar dt, const DiphasicStepData& stepData ,
	Phase& phase ) const
{
	typedef Eigen::SparseMatrix< Scalar > SM ;

	(void) dt ;

	// Compute rhs of momentum conservation -- gravity + u(t)
	DynVec rhs ;
	{
		DynVec forces = stepData.forms.externalForces ;

		// Inertia
#ifndef  FULL_FEM
		forces += stepData.forms.linearMomentum  ;
#endif

		rhs = stepData.fullGridProj.vel * forces ;
	}

	// Solve unconstrained momentum equation
	const Index m  = stepData.forms.A.rows() ;
	const Index ma = stepData.forms.R.rows() ;
	const Index n  = stepData.forms.B.rows() ;

	DynVec x ( m+ma+n ) ;
	x.head( m ).setZero() ; // = stepData.forms.M_lumped_inv * rhs ;
	x.segment( m, ma ).setZero() ;
	x.tail(n).setZero();

	DynVec l ( m+ma+n ) ;
	l.head(m) = rhs ;
	l.segment( m, ma ).setZero() ;
	l.tail(n).setZero();

	SM M ;
	makePenalizedEigenStokesMatrix( stepData, M );
	Eigen::SimplicialLDLT< SM > M_fac ( M ) ;
	Log::Verbose() << "LU factorization : " << M_fac.info() << std::endl ;

	x = M_fac.solve( l ) ;

	DynVec wvh ( WD * stepData.nPrimalNodes() ) ;
	wvh.setZero() ;

	// Output
	PrimalVectorField w( phase.velocity.shape() ) ;
	stepData.primalNodes.var2field( wvh, w ) ;

	// U_1
	phase.velocity.flatten() = x.head(m) + w.flatten() ;

	// U_2
	PrimalScalarField ratio( phase.fraction.shape() ) ;
	ratio.flatten() = phase.fraction.flatten().array() / (1. - phase.fraction.flatten().array() ) ;
	w.multiply_by( ratio ) ;

	phase.geo_proj.flatten() = x.head(m) - (config.alpha()+1)*w.flatten() ;
}


} // d6
