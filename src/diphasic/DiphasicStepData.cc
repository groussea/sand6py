#include "DiphasicStepData.hh"

#include "FluidPhase.hh"

#include "mono/Phase.hh"
#include "mono/PhaseStepData.hh"

#include "simu/FormBuilder.impl.hh"
#include "simu/DynParticles.hh"

#include "geo/BoundaryInfo.hh"

#include "utils/Config.hh"
#include "utils/Log.hh"

#include <bogus/Core/Block.impl.hpp>
#include <bogus/Core/Utils/Timer.hpp>

//#define KINEMATIC_VISC
//#define U_MOMENTUM


// a C' p  = a wh / Stk
// a C  wh + Bu = 0

// C'Stk p  = Stk (wh/Stk)
// C Stk wh/Stk + Bu = 0

// C'sStk p  = sStk/a (a wh/Stk) = 1/a (a wh/sStk)
// C sStk (a wh/sStk) + Bu = 0

//u1 = u + wh = u + sStk/a (1 wh/sStk)

namespace d6 {

const Scalar DiphasicStepData::s_maxPhi = 0.999 ;

void DiphasicStepData::computeProjectors(const Config&config,
									  const PrimalShape& pShape, FullProjectors& mats )
{
	//Full grid projectors

	const Index m  = pShape.nDOF() ;

	mats.vel.setRows( m );
	mats.vel.setIdentity() ;

	mats.pressure.setRows( m );
	mats.pressure.setIdentity() ;

	StrBoundaryMapper bdMapper ( config.boundary ) ;

	MatS stressProj ;

	for( auto cellIt = pShape.mesh().cellBegin() ;
		 cellIt != pShape.mesh().cellEnd() ; ++cellIt )
	{
		if( ! pShape.mesh().onBoundary( *cellIt ) ) continue ;

		typename PrimalShape::Location ploc ;
		typename PrimalShape::NodeList pnodes ;
		ploc.cell = *cellIt ;
		pShape.list_nodes( ploc, pnodes ) ;

		for( unsigned k = 0 ; k < PrimalShape::NI ; ++k ) {
			const Index i = pnodes[k] ;

			BoundaryInfo info ;
			pShape.locate_dof( ploc, k ) ;
			pShape.mesh().boundaryInfo( ploc, bdMapper, info ) ;
			info.velProj( mats.vel.block( i ) ) ;

			if( !config.weakStressBC ) {
				info.stressProj( stressProj ) ;
				mats.pressure.block( i ) = stressProj.block<1,1>(0,0) ;
			}
		}

	}

}

void DiphasicStepData::assembleMatrices(
		const Particles &particles,
		const Config &config, const Scalar dt, const DualShape &dShape,
		const FluidPhase& fluid, const Phase& phase,
		const PrimalScalarField &intPhi, const PrimalVectorField& intPhiVel )
{
	const Scalar mass_regul = 1.e-8 ;
	const Scalar sStk = 1./std::sqrt( config.fluidFriction ) ;

	bogus::Timer timer ;

	const PrimalShape& pShape = intPhi.shape() ;

	typedef typename PrimalShape::Interpolation P_Itp ;
	typedef typename PrimalShape::Derivatives   P_Dcdx ;
	typedef typename DualShape::Interpolation   D_Itp ;
	typedef typename DualShape::Derivatives     D_Dcdx ;


	// Projectors
	   PhaseStepData::computeProjectors( config, pShape, dShape, intPhi, primalNodes, dualNodes, 0, activeProj ) ;
	DiphasicStepData::computeProjectors( config, pShape, fullGridProj ) ;

	const Index m  = pShape.nDOF() ;


	std::vector<Index> fullIndices( m );
	std::iota( fullIndices.begin(), fullIndices.end(), 0);

	// I - Full grid forms (A,B)
	{


		typedef FormBuilder< PrimalShape, PrimalShape > Builder ;
		Builder builder( pShape, pShape ) ;
		builder.reset( m );

		builder.addToIndex( fullIndices, fullIndices );
		builder.makeCompressed();

		forms.A.clear();
		forms.A.setRows( m );
		forms.A.setCols( m );
		forms.A.cloneIndex( builder.index() ) ;
		forms.A.setBlocksToZero() ;

		forms.B.clear();
		forms.B.setRows( m );
		forms.B.setCols( m );
		forms.B.cloneIndex( builder.index() ) ;
		forms.B.setBlocksToZero() ;

		builder.integrate_qp( [&]( Scalar w, const Vec&, const P_Itp& l_itp, const P_Dcdx& l_dc_dx, const P_Itp& r_itp, const P_Dcdx& r_dc_dx )
		{
			Builder:: addDuDv( forms.A, w*2*config.viscosity, l_itp, l_dc_dx, r_itp, r_dc_dx, fullIndices, fullIndices ) ;
			Builder:: addDpV ( forms.B, -w, l_itp, l_dc_dx, r_itp, fullIndices, fullIndices ) ;
		}
		);
		Log::Debug() << "A Integrate grid: " << timer.elapsed() << std::endl ;


		// Linear forms
		{

			// Integrated quantities
			PrimalVectorField linearMomentum( pShape ) ;
			PrimalScalarField beta( pShape ) ;

			linearMomentum.flatten() = (1+config.alpha()) * intPhiVel.flatten() ;
			beta.flatten()           = config.alpha() * intPhi.flatten() ;

#ifdef U_MOMENTUM
			linearMomentum.set_zero() ;
#endif

			typename PrimalShape::Interpolation itp ;
			typename PrimalShape::Location loc ;

//			builder.integrate_node( pShape.mesh().cellBegin(), pShape.mesh().cellEnd(),
//									[&]( Scalar w, const Vec&pos, const P_Itp& itp, const P_Itp& ) {
			builder.integrate_qp( [&]( Scalar w, const Vec&pos, const P_Itp& itp, const P_Dcdx&, const P_Itp&, const P_Dcdx&) {

				const Scalar phi = std::min( s_maxPhi, phase.fraction( pos ) ) ;
#ifdef U_MOMENTUM
				const Vec u = fluid.mavg_vel( pos ) ;
				Vec pos_prev = pShape.mesh().clamp_point( pos - dt*u ) ;
				const Vec u_adv = fluid.mavg_vel( pos_prev ) ;
#else
				const Vec u2 = fluid.velocity( pos ) ;
				Vec pos_prev = pShape.mesh().clamp_point( pos - dt*u2 ) ;
				const Vec u2_adv = fluid.velocity( pos_prev ) ;
#endif
				for( Index k = 0 ; k < itp.nodes.size() ; ++k ) {
					beta[ itp.nodes[k] ]           += w * itp.coeffs[k] ;
#ifdef U_MOMENTUM
					linearMomentum[ itp.nodes[k] ] += w * itp.coeffs[k] * u_adv ;
#else
					linearMomentum[ itp.nodes[k] ] += w * itp.coeffs[k] * (1 - phi) * u2_adv ;
#endif
				}
			} ) ;

#ifdef U_MOMENTUM
			builder.integrate_particle( particles, [&]( Index i, Scalar w, const P_Itp& itp, const P_Dcdx& , const P_Itp& , const P_Dcdx& )
			{
				const Vec& pos = particles.centers().col(i) ;
				const Vec u = fluid.mavg_vel( pos ) ;
				Vec pos_prev = pShape.mesh().clamp_point( pos - dt*u ) ;
				const Vec u_adv = fluid.mavg_vel( pos_prev ) ;
				for( Index k = 0 ; k < itp.nodes.size() ; ++k ) {
					linearMomentum[ itp.nodes[k] ] += w * itp.coeffs[k] * config.alpha() * u_adv ;
				}
			} ) ;
#endif

			forms.linearMomentum = linearMomentum.flatten() * config.fluidVolMass / dt ;

			PrimalVectorField gravity ( pShape ) ;
			gravity.set_constant( config.gravity );
			gravity.multiply_by( beta ) ;
			gravity.flatten() *= config.fluidVolMass  ;

			forms.externalForces = gravity.flatten() ;



			// Lumped mass matrix
			{
				forms.M_lumped.setRows( m );
				forms.M_lumped.setIdentity() ;
				forms.M_lumped_inv.setRows( m );
				forms.M_lumped_inv.setIdentity() ;
	#pragma omp parallel for
				for( Index i = 0 ; i < m ; ++i ) {
					forms.M_lumped.block( i ) *= beta[ i ] * config.fluidVolMass /  dt ;
				}
			}
		}

#ifndef FULL_FEM
		forms.A += forms.M_lumped ;
#endif

		// Projections
		const typename FormMat<WD,WD>::SymType IP = fullGridProj.vel.Identity() - fullGridProj.vel ;
		forms.A = fullGridProj.vel * forms.A * fullGridProj.vel  + IP ;

#pragma omp parallel for
		for( Index i = 0 ; i < m ; ++i ) {
			const Scalar mass = forms.M_lumped.block(i).trace() / WD ;
			forms.M_lumped         .block(i) = fullGridProj.vel.block(i) * mass
					+ Mat::Identity() - fullGridProj.vel.block(i) ;
			forms.M_lumped_inv     .block(i) = fullGridProj.vel.block(i) * 1./(mass + mass_regul )
					+ Mat::Identity() - fullGridProj.vel.block(i) ;
		}

	}


	// II Active forms -- R, C

	{
		const Index ma = primalNodes.count() ;

		typedef FormBuilder< PrimalShape, PrimalShape > Builder ;
		Builder builder( pShape, pShape ) ;
		builder.reset( m );

		builder.addToIndex( fullIndices, primalNodes.indices );
		builder.makeCompressed();

		forms.R.clear();
		forms.R.setRows( ma );
		forms.R.setCols( ma );
		forms.R.setIdentity() ;

		forms.C.clear();
		forms.C.setRows( m );
		forms.C.setCols( ma );
		forms.C.cloneIndex( builder.index() ) ;
		forms.C.setBlocksToZero() ;

		DynVec Rcoeffs ( ma ) ;
		Rcoeffs.setZero() ;


		builder.integrate_particle( particles, [&]( Index i, Scalar w, const P_Itp& l_itp, const P_Dcdx& l_dc_dx, const P_Itp& r_itp, const P_Dcdx& r_dc_dx)
		{
			Builder:: addDpV ( forms.C, w*sStk, l_itp, l_dc_dx, r_itp, fullIndices, primalNodes.indices ) ;

#ifdef KINEMATIC_VISC
			Builder:: addDuDv( forms.A, w*config.alpha()*2*config.viscosity, l_itp, l_dc_dx, r_itp, r_dc_dx, fullIndices, fullIndices ) ;
#else
			(void) r_dc_dx ;
#endif

			// TODO: test w/ other funcs
			const Vec& pos = particles.centers().col(i) ;
			const Scalar phi = std::min( s_maxPhi, phase.fraction( pos ) ) ;
			const Scalar vR = (1 + config.alpha() * phi )/(1-phi) ;
//			const Scalar vR = 1;

			for( Index k = 0 ; k < l_itp.nodes.size() ; ++k ) {
				const Index idx = primalNodes.indices[ l_itp.nodes[k]] ;
				Rcoeffs[ idx ] += w * l_itp.coeffs[k] / config.alpha() * vR ;
			}

		}
		);
		Log::Debug() << "C Integrate particle: " << timer.elapsed() << std::endl ;

		{
#pragma omp parallel for
			for( Index i = 0 ; i < ma ; ++i ) {
				forms.R.block(i) = activeProj.vel.block(i) * Rcoeffs[i]
						+ Mat::Identity() - activeProj.vel.block(i) ;
			}
		}

	}


	const Index n = dualNodes.count() ;
	forms.S.compute( dShape, dualNodes, n );

	// III. G
	{
		typedef FormBuilder< DualShape, PrimalShape > Builder ;
		Builder builder( dShape, pShape ) ;
		builder.reset( n );

		builder.addToIndex( dualNodes.indices, fullIndices );
		builder.makeCompressed();

		forms.G.clear();
		forms.G.setRows( n );
		forms.G.setCols( m );
		forms.G.cloneIndex( builder.index() ) ;
		forms.G.setBlocksToZero() ;

		builder.integrate_particle( particles, [&]( Index, Scalar w, const D_Itp& l_itp, const D_Dcdx& , const P_Itp& r_itp, const P_Dcdx& r_dc_dx )
		{
			Builder::addTauDu( forms.G, w, l_itp, r_itp, r_dc_dx, dualNodes.indices, fullIndices ) ;
		} ) ;
	}

	// IV. H
	{
		const Index ma = primalNodes.count() ;

		typedef FormBuilder< DualShape, PrimalShape > Builder ;
		Builder builder( dShape, pShape ) ;
		builder.reset( n );

		builder.addToIndex( dualNodes.indices, primalNodes.indices );
		builder.makeCompressed();

		forms.H.clear();
		forms.H.setRows( n );
		forms.H.setCols( ma );
		forms.H.cloneIndex( builder.index() ) ;
		forms.H.setBlocksToZero() ;

		builder.integrate_particle( particles, [&]( Index, Scalar w, const D_Itp& l_itp, const D_Dcdx& , const P_Itp& r_itp, const P_Dcdx& r_dc_dx )
		{
			Builder::addTauDu( forms.H, w*sStk/config.alpha(), l_itp, r_itp, r_dc_dx, dualNodes.indices, primalNodes.indices ) ;
		} ) ;
	}

}




void DiphasicStepData::compute(const DynParticles& particles,
		const Config &config, const Scalar dt, const FluidPhase& fluid, Phase &phase )
{

	const PrimalShape& pShape = phase.velocity.shape() ;
	const   DualShape& dShape = phase.stresses.shape() ;

	// Transfer particles quantities to grid
	PrimalScalarField intPhiPrimal  ( pShape ) ;
	PrimalVectorField intPhiVel     ( pShape ) ;

	// Not used, but required by particle interface
	DualScalarField intPhiDual    ( dShape ) ;
	DualScalarField intPhiInertia ( dShape ) ;
	DualScalarField intPhiCohesion( dShape ) ;
	DualTensorField intPhiOrient  ( dShape ) ;

	std::vector< bool > activeCells ;

#if defined(FULL_FEM)
	pShape.compute_tpz_mass   ( intPhiPrimal.flatten() ) ;
	dShape.compute_lumped_mass( intPhiDual  .flatten() ) ;
	intPhiVel.set_zero() ;
	intPhiInertia.set_zero() ;
	intPhiCohesion.set_zero() ;
	intPhiOrient.set_zero() ;
	activeCells.assign( pShape.mesh().nCells(), true ) ;
#else
	particles.integratePrimal( activeCells, intPhiPrimal, intPhiVel ) ;
	particles.integrateDual( intPhiDual, intPhiInertia, intPhiOrient, intPhiCohesion ) ;
#endif

	// Compute phi and grad_phi (for visualization purposes )
	PhaseStepData::computePhiAndGradPhi( intPhiPrimal, phase.fraction, phase.grad_phi ) ;
	Log::Debug() << "MAX phi " << phase.fraction.max_abs() << std::endl ;

	phase.velocity = intPhiVel  ;
	phase.velocity.divide_by_positive( intPhiPrimal ) ;

	// Active nodes
	PhaseStepData::computeActiveNodes( activeCells, pShape, dShape, primalNodes, dualNodes ) ;
	Log::Verbose() << "Active nodes: " << nPrimalNodes() << " / " << pShape.nDOF() << std::endl;
	Log::Verbose() << "  Dual nodes: " <<   nDualNodes() << " / " << dShape.nDOF() << std::endl;

	// Bilinear forms matrices
	assembleMatrices( particles.geo(), config, dt, dShape, fluid, phase, intPhiPrimal, intPhiVel );

	dualNodes  .field2var( intPhiDual, forms.fraction ) ;

	// Cohesion inertia
	intPhiCohesion.divide_by_positive( intPhiDual ) ;
	intPhiInertia .divide_by_positive( intPhiDual ) ;
	dualNodes.field2var( intPhiCohesion, cohesion ) ;
	dualNodes.field2var( intPhiInertia , inertia  ) ;

	// Volumes
	{
		DualScalarField volumes ( dShape ) ;
		dShape.compute_lumped_mass( volumes.flatten() );
		dualNodes.field2var( volumes, forms.volumes ) ;
	}
}

} // d6
