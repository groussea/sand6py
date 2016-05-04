#include "DiphasicSolver.hh"

#include "DiphasicStepData.hh"

#include "FluidPhase.hh"
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
		Bproj = stepData.fullGridProj.pressure * ( stepData.forms.B * stepData.fullGridProj.vel ) ;
		bogus::convert( Bproj, B ) ;

		typename FormMat<1,WD>::Type Cproj ;
		Cproj = stepData.fullGridProj.pressure * ( stepData.forms.C * stepData.activeProj.vel ) ;
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

void DiphasicSolver::step(const Config &config, const Scalar dt, Phase &phase, FluidPhase &fluid ) const
{

	DiphasicStepData stepData ;
	stepData.compute( m_particles, config, dt, fluid, phase );

	solve( config, dt, stepData, phase, fluid ) ;
}


void DiphasicSolver::solve(
	const Config& config, const Scalar dt, const DiphasicStepData& stepData ,
	Phase& phase, FluidPhase &fluid ) const
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
	x.head(m) = phase.velocity.flatten() ; // = stepData.forms.M_lumped_inv * rhs ;
	x.tail(n) = fluid.pressure.flatten() ;

	auto w = x.segment( m, ma ) ;
	stepData.primalNodes.field2var( fluid.velocity, w ) ;

	DynVec l ( m+ma+n ) ;
	l.head(m) = rhs ;
	l.segment( m, ma ).setZero() ;
	l.tail(n).setZero();

	SM M ;
	makePenalizedEigenStokesMatrix( stepData, M );
	Eigen::SimplicialLDLT< SM > M_fac ( M ) ;
	Log::Verbose() << "LU factorization : " << M_fac.info() << std::endl ;

	x = M_fac.solve( l ) ;

	// Output
	fluid.pressure.flatten() = x.tail(n) ;

	// U_1 = u + wh
	stepData.primalNodes.var2field( w, fluid.velocity ) ;
	phase.velocity.flatten() = x.head(m) + fluid.velocity.flatten() ;

	// U_2 = u - (alpha+1) phi/(1-phi) wh
	PrimalScalarField ratio( phase.fraction.shape() ) ;
	ratio.flatten() = phase.fraction.flatten().array() / (1. - phase.fraction.flatten().array().min(config.phiMax) ) ;
	fluid.velocity.multiply_by( ratio ) ;

	fluid.velocity.flatten() = x.head(m) - (config.alpha()+1)*fluid.velocity.flatten() ;

	std::cout << "U " << x.head(m).lpNorm< Eigen::Infinity >() << std::endl ;
	std::cout << "W " << ( fluid.velocity.flatten() - phase.velocity.flatten()).lpNorm< Eigen::Infinity >() << std::endl ;
}


} // d6
