#include "DiphasicStepData.hh"

#include "FluidPhase.hh"

#include "mono/Phase.hh"
#include "mono/PhaseStepData.hh"

#include "simu/FormBuilder.impl.hh"
#include "simu/DynParticles.hh"

#include "geo/geo.fwd.hh"
#include "geo/BoundaryInfo.hh"
#include "geo/FieldBase.impl.hh"

#include "utils/Config.hh"
#include "utils/Log.hh"

#include <bogus/Core/Block.impl.hpp>
#include <bogus/Core/Utils/Timer.hpp>

#include <functional>
#include <numeric>

#define  W_VISC
#define WU_VISC
#define W_MOMENTUM

#define U2_LUMPING (.0)
#define  W_LUMPING (1.)

#if D6_MESH_IMPL != D6_MESH_TET_GRID && D6_DIM == 2
#define INTEGRATE_LOBATTO
#endif

#define EI_VISC( a )( 1 + 2.5*(a) )

// a C' p  = a wh / Stk
// a C  wh + Bu = 0

// C'Stk p  = Stk (wh/Stk)
// C Stk wh/Stk + Bu = 0

// C'sStk p  = sStk/a (a wh/Stk) = 1/a (a wh/sStk)
// C sStk (a wh/sStk) + Bu = 0

//u1 = u + wh = u + sStk/a (1 wh/sStk)

// 1/(dt) (w) = a w / Stk
// wh = a sStk


// Stk (1-phi)/beta w = wh sStk /a
// wh = a * sStk (1-phi)/beta
// w  =  beta/a(1-phi)  wh/sStk

// a phi C' p = phi Stk/dt dw
// sStk C' p = sStk Stk phi/(a dt) dw
//    = phi beta/(1-phi)  Stk /(dt a^2) dwh
//    = phi beta/a(1-phi)  Stk/(dt a) wh  - \phi Stk sStk/(a dt) wk
// R += phi beta/a(1-phi) (Stk /dt a)


// a Stk/Re \phi Div( D(a pi /beta w))
// = a sStk/Re \phi Div( D(phi wh))
// (* sStk/a) = Stk/Re \phi Div( D(phi wh))

// 1/(dt) (w - wk) = a w / Stk
// 1/(dt) (w - wk) = a w / Stk

// 1/(dt) (wh - (1-phi)/beta wk) = a wh / Stk

// phi/(a sStk dt) (dw) = a wh / Stk
// (1-phi)/beta 1/(dt) (w - wk) = a wh / Stk

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
	const Scalar w_regul = 1.e-2 / config.alpha() ; // Small w.r.t phi max (1./alpha, phi)
	const Scalar sStk = std::sqrt( config.Stokes() ) ;

	//Note: config.viscosity = 1./(Rey (alpha+1))

	bogus::Timer timer ;

	const PrimalShape& pShape = intPhi.shape() ;

	typedef typename PrimalShape::Interpolation P_Itp ;
	typedef typename PrimalShape::Derivatives   P_Dcdx ;
	typedef typename DualShape::Interpolation   D_Itp ;
	typedef typename DualShape::Derivatives     D_Dcdx ;


	// Projectors
	   PhaseStepData::computeProjectors( config, pShape, dShape, intPhi, primalNodes, dualNodes, 0, activeProj ) ;
	DiphasicStepData::computeProjectors( config, pShape, fullGridProj ) ;


	// Characteristics functor
	// TODO subsample
	const auto char_func = [dt,&pShape]( const PrimalVectorField &velocity, const Vec &pos ) {
		const Vec u2 = velocity( pos ) ;
		const Vec pos_prev = pShape.mesh().clamp_point( pos - dt*u2 ) ;
		const Vec u2_prev  = velocity( pos_prev ) ;
		const Vec trial = pos - .5*dt*(u2+u2_prev) ;
		const Vec clamped =  pShape.mesh().clamp_point( trial ) ;
		Vec u2_adv = velocity( clamped ) ;
		const Scalar projn2 = (clamped - trial).squaredNorm() ;
		if( projn2 > 1.e-6 ) {
			u2_adv -= u2_adv.dot(clamped-trial)/projn2*(clamped-trial) ;
		}
		return u2_adv ;
	} ;

	const auto u2_adv_func = std::bind( char_func, fluid.velocity, std::placeholders::_1 ) ;

	PrimalVectorField u2_adv_field ( pShape ) ;
	u2_adv_field.set_free( u2_adv_func ) ;


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

		// 1/((alpha+1) Rey) D(u):D(v) + (beta)/dt
		forms.A.clear();
		forms.A.setRows( m );
		forms.A.setCols( m );
		forms.A.cloneIndex( builder.index() ) ;
		forms.A.setBlocksToZero() ;

		// 1 q div v
		forms.B.clear();
		forms.B.setRows( m );
		forms.B.setCols( m );
		forms.B.cloneIndex( builder.index() ) ;
		forms.B.setBlocksToZero() ;

		// D(u):tau (for computing strain rate)
		forms.D.clear();
		forms.D.setRows( m );
		forms.D.setCols( m );
		forms.D.cloneIndex( builder.index() ) ;
		forms.D.setBlocksToZero() ;

		builder.integrate_qp( [&]( Scalar w, const Vec& pos, const P_Itp& l_itp, const P_Dcdx& l_dc_dx, const P_Itp& r_itp, const P_Dcdx& r_dc_dx )
		{
			const Scalar visc = 2 * config.viscosity * EI_VISC( phase.fraction( pos ) ) ;
			Builder:: addDuDv ( forms.A, w*visc, l_itp, l_dc_dx, r_itp, r_dc_dx, fullIndices, fullIndices ) ;
			Builder:: addQDivu( forms.B, w, l_itp, r_itp, r_dc_dx, fullIndices, fullIndices ) ;
			Builder:: addTauDu( forms.D, w, l_itp, r_itp, r_dc_dx, fullIndices, fullIndices ) ;
		}
		);
		Log::Debug() << "A Integrate grid: " << timer.elapsed() << std::endl ;


		// Linear forms
		{

			// Integrated quantities
			PrimalVectorField linearMomentum( pShape ) ;
			PrimalScalarField beta( pShape ) ;

			linearMomentum.flatten() = (1+config.alpha()) * intPhiVel.flatten() ;
			beta.flatten()           = (1+config.alpha()) * intPhi.flatten() ;

#ifdef U_MOMENTUM
			linearMomentum.set_zero() ;
#endif

			typename PrimalShape::Interpolation itp ;
			typename PrimalShape::Location loc ;

#ifdef INTEGRATE_LOBATTO
			builder.integrate_node( pShape.mesh().cellBegin(), pShape.mesh().cellEnd(),
			                        [&]( Scalar w, const Vec&pos, const P_Itp& itp, const P_Itp& ) {

			w *= PrimalShape::NI * 1./36  ;
#else
			builder.integrate_qp( [&]( Scalar w, const Vec&pos, const P_Itp& itp, const P_Dcdx&, const P_Itp&, const P_Dcdx&) {
#endif //LOBATTO

				const Scalar phi = std::min(1., phase.fraction( pos ) ) ;
#ifdef U_MOMENTUM
				const Vec u_adv = char_func( fluid.mavg_vel, pos ) ;
#else
				const Vec u2_adv = u2_adv_func( pos ) ;
#endif // U_MOMENTUM
				for( Index k = 0 ; k < itp.nodes.size() ; ++k ) {
					const Scalar weight =  w * itp.coeffs[k] * (1 - phi) ;
					beta[ itp.nodes[k] ]           += weight ;
#ifdef U_MOMENTUM
					linearMomentum[ itp.nodes[k] ] += w * itp.coeffs[k] * u_adv ;
#else
					linearMomentum[ itp.nodes[k] ] += (1-U2_LUMPING) * weight * u2_adv ;
					linearMomentum[ itp.nodes[k] ] += (  U2_LUMPING) * weight * u2_adv_field[ itp.nodes[k] ] ;

					const Scalar base = (1-U2_LUMPING) * weight * config.fluidVolMass / dt ;
					for( Index j = 0 ; j < itp.nodes.size() ; ++j ) {
						forms.A.block( itp.nodes[k], itp.nodes[j] ).diagonal() += Vec::Constant(base * itp.coeffs[j]) ;
					}
#endif // U_MOMENTUM

				}

			} ) ;



#ifdef INTEGRATE_LOBATTO

			for( auto cellIt = pShape.mesh().cellBegin() ; cellIt != pShape.mesh().cellEnd() ; ++cellIt )
			{
				PrimalMesh::CellGeo geo ;
				pShape.mesh().get_geo( *cellIt, geo );
				PrimalMesh::Location loc ;
				loc.cell = *cellIt ;

				DynMatW clist( WD, 5 ) ;
				clist <<  0,.5,.5,.5,1.,
				         .5,0 ,.5, 1,0 ;
				DynArr  wlist(5) ;
				wlist <<  1./9,1./9,4./9,1./9, 1./9 ;

				for(unsigned k = 0 ; k < clist.cols() ; ++k ) {
				const Scalar w = geo.volume() * wlist[k] ;

				loc.coords = clist.col(k) ; //FIXME TET
				PrimalShape::Interpolation itp ;
				pShape.interpolate( loc, itp );
				const Vec& pos = pShape.mesh().pos( loc ) ;

				const Scalar phi = std::min(1., phase.fraction( pos ) ) ;
				const Vec u2_adv = u2_adv_func( pos ) ;

				for( Index k = 0 ; k < itp.nodes.size() ; ++k ) {
					const Scalar weight =  w * itp.coeffs[k] * (1 - phi) ;
					beta[ itp.nodes[k] ]           += weight ;
					linearMomentum[ itp.nodes[k] ] += (1-U2_LUMPING) * weight * u2_adv ;
					linearMomentum[ itp.nodes[k] ] += (  U2_LUMPING) * weight * u2_adv_field[ itp.nodes[k] ] ;

					const Scalar base = (1-U2_LUMPING) * weight * config.fluidVolMass / dt ;
					for( Index j = 0 ; j < itp.nodes.size() ; ++j ) {
						forms.A.block( itp.nodes[k], itp.nodes[j] ).diagonal() += Vec::Constant(base * itp.coeffs[j]) ;
					}

				}
				}

			}
#endif
			// At this point A = Avisc + (1-U2_LUMPING) rho_f/dt (1-\phi)

			forms.linearMomentum = linearMomentum.flatten() * config.fluidVolMass / dt ;

			PrimalVectorField gravity ( pShape ) ;
			gravity.set_constant( config.gravity );
			gravity.multiply_by( beta ) ;
			gravity.flatten() *= config.fluidVolMass  ;

			forms.externalForces = gravity.flatten() ;



			// Lumped mass matrix
			{
    #pragma omp parallel for
				for( Index i = 0 ; i < m ; ++i ) {
					forms.A.block( i,i ).diagonal() += Vec::Constant( (1-U2_LUMPING)*( (config.alpha() + 1) * intPhi[i] * config.fluidVolMass /  dt + mass_regul ) );
				}
				// At this point A = Avisc + (1-U2_LUMPING) rho_f/dt \beta

				forms.M_lumped.setRows( m );
				forms.M_lumped.setIdentity() ;
				forms.M_lumped_inv.setRows( m );
				forms.M_lumped_inv.setIdentity() ;
    #pragma omp parallel for
				for( Index i = 0 ; i < m ; ++i ) {
					forms.M_lumped.block( i ).diagonal() =  Vec::Constant( beta[ i ] * config.fluidVolMass /  dt );
				}
				forms.A += U2_LUMPING * forms.M_lumped ;
			}
		}

		// Wind, increasing with altitude
		PrimalVectorField bdVel( pShape ) ;
		bdVel.set_free( [&config] (const Vec &x){
//			return config.windSpeed * ( 2 - x[WD-1]/config.box[WD-1])*(x[WD-1]/config.box[WD-1]) ;
			return config.windSpeed * (1. -  std::exp( - 10 * x[WD-1]/config.box[WD-1] ) ) ;
		} ) ;

		// Projections
		const typename FormMat<WD,WD>::SymType IP = fullGridProj.vel.Identity() - fullGridProj.vel ;

		dirichletVel = IP * bdVel.flatten() ;
		forms.dirichletTerm = dirichletVel + fullGridProj.vel * forms.A * dirichletVel;

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

	{
		const Index ma = primalNodes.count() ;

		typedef FormBuilder< PrimalShape, PrimalShape > Builder ;
		Builder builder( pShape, pShape ) ;
		builder.reset( ma );

		builder.addToIndex( primalNodes.indices, primalNodes.indices );
		builder.makeCompressed();

		// 1/(Rey (alpha+1)) (\alpha^2) D(\phi w) D(\phi z)
		forms.R_visc.clear();
		forms.R_visc.setRows( ma );
		forms.R_visc.setCols( ma );
		forms.R_visc.finalize();
		forms.R_visc.cloneIndex( builder.index() ) ;
		forms.R_visc.setBlocksToZero() ;

#ifdef W_VISC

		builder.integrate_cell<form::Left>( primalNodes.cells.begin(), primalNodes.cells.end(),
		                        [&]( Scalar w, const Vec& pos, const P_Itp& l_itp, const P_Dcdx& l_dc_dx, const P_Itp& r_itp, const P_Dcdx& )
		    {
			    // u = sum_k ( c_k uK  )
			    // du_dx  = sum_k ( dcdx_k  uK  )
			    // phiu   = sum_k ( c_k phiK uK  )
			    // d_phiu = sum_k ( (dcdx_k phiK) uK  )
			    P_Dcdx phidc_dx ;
				for( Index k = 0 ; k < l_itp.nodes.size() ; ++k ) {
					phidc_dx.row(k) = l_dc_dx.row(k) * phase.fraction[ l_itp.nodes[k] ] ;
				}
				const Scalar visc = 2 * config.viscosity * EI_VISC( phase.fraction( pos ) ) ;
				const Scalar weight = w * visc * config.Stokes() * config.alpha() * config.alpha() ;
				Builder:: addDuDv( forms.R_visc, weight, l_itp, phidc_dx, r_itp, phidc_dx, primalNodes.indices, primalNodes.indices ) ;
		    }
		);
#endif

	}

	// II Active forms -- R, C

	{
		const Index ma = primalNodes.count() ;

		typedef FormBuilder< PrimalShape, PrimalShape > Builder ;
		Builder builder( pShape, pShape ) ;
		builder.reset( m );

		builder.addToIndex( fullIndices, primalNodes.indices );
		builder.makeCompressed();

		// sStk \alpha / (\Rey (alpha +1)) D(phi w)D(v)
		forms.F.clear();
		forms.F.setRows( m );
		forms.F.setCols( ma );

#ifdef WU_VISC
		forms.F.cloneIndex( builder.index() ) ;
		forms.F.setBlocksToZero() ;

		builder.integrate_cell<form::Left>( primalNodes.cells.begin(), primalNodes.cells.end(),
		                        [&]( Scalar w, const Vec&pos, const P_Itp& l_itp, const P_Dcdx& l_dc_dx, const P_Itp& r_itp, const P_Dcdx& )
		    {
			    // u = sum_k ( c_k uK  )
			    // du_dx  = sum_k ( dcdx_k  uK  )
			    // phiu   = sum_k ( c_k phiK uK  )
			    // d_phiu = sum_k ( (dcdx_k phiK) uK  )
			    P_Dcdx phidc_dx ;
				for( Index k = 0 ; k < l_itp.nodes.size() ; ++k ) {
					phidc_dx.row(k) = l_dc_dx.row(k) * phase.fraction[ l_itp.nodes[k] ] ;
				}

				const Scalar visc = 2 * config.viscosity * EI_VISC( phase.fraction( pos ) ) ;

				const Scalar weight = - w * visc * sStk * config.alpha() ;
				Builder:: addDuDv( forms.F, weight, l_itp, l_dc_dx, r_itp, phidc_dx, fullIndices, primalNodes.indices ) ;
		    }
		);
#endif

		// phi beta / (1-\phi)  (St/dt + beta \Xi )
		forms.R.clear();
		forms.R.setRows( ma );
		forms.R.setCols( ma );
		forms.R.setIdentity() ;

		// sStk * alpha \phi grad_p . v
		forms.C.clear();
		forms.C.setRows( m );
		forms.C.setCols( ma );
		forms.C.cloneIndex( builder.index() ) ;
		forms.C.setBlocksToZero() ;

		DynVec Rcoeffs ( ma ) ;
		Rcoeffs.setZero() ;

		PrimalVectorField fluctuMomentum( pShape ) ;
		fluctuMomentum.set_zero() ;

		const Scalar St_dt  = config.Stokes() / dt ;
		Log::Debug() << "Ratio inertia flutuation " << St_dt * std::pow(1 - config.phiMax, config.RZExponent) / (config.alpha()+1) << std::endl ;

		builder.integrate_particle( particles, [&]( Index i, Scalar w, const P_Itp& l_itp, const P_Dcdx& l_dc_dx, const P_Itp& r_itp, const P_Dcdx& r_dc_dx)
		{
			Builder:: addDpV ( forms.C, w*sStk*config.alpha(), l_itp, l_dc_dx, r_itp, fullIndices, primalNodes.indices ) ;

			// TODO: test w/ other funcs
			const Vec& pos = particles.centers().col(i) ;
			const Scalar phi = phase.fraction( pos ) ;
			const Scalar mphi = (1-std::min( config.phiMax, phi )) ;
			const Scalar beta = (config.alpha()*phi + 1) ;
			const Scalar beta_mphi = config.volMass * beta/mphi ;
			Vec vW = Vec::Zero() ;

			const Scalar RZ_54 = std::pow(mphi, -config.RZExponent ) ; // Ridcharson-Zaki closure
			Scalar vR = beta_mphi * ( config.alpha() * beta )/(config.alpha() + 1) * RZ_54;
#ifdef W_MOMENTUM
			if( phi < s_maxPhi ) {
				const Vec& u1_adv = particles.velocities().col(i) ;
				const Vec& u2_adv = u2_adv_func( pos ) ;

				vW  = St_dt/sStk * ( u1_adv - u2_adv ) ;
			}

			vR += beta_mphi * St_dt ;
#endif

			for( Index k = 0 ; k < l_itp.nodes.size() ; ++k ) {
				const Index idx = primalNodes.indices[ l_itp.nodes[k]] ;
				fluctuMomentum[ l_itp.nodes[k] ] += w * l_itp.coeffs[k] * vW;

				const Scalar base = w * l_itp.coeffs[k] * vR ;
				Rcoeffs[ idx ] += base;

				for( Index j = 0 ; j < l_itp.nodes.size() ; ++j ) {
					forms.R_visc.block( idx, primalNodes.indices[l_itp.nodes[j]] ).diagonal() += Vec::Constant( (1 - W_LUMPING) * base * l_itp.coeffs[j]) ;
				}
			}

		}
		);
		Log::Debug() << "C Integrate particle: " << timer.elapsed() << std::endl ;

		primalNodes.field2var( fluctuMomentum, forms.fluctuMomentum ) ;

		{
#pragma omp parallel for
			for( Index i = 0 ; i < ma ; ++i ) {
				forms.R.block(i) = activeProj.vel.block(i) * ( Rcoeffs[i] + w_regul )
				        + Mat::Identity() - activeProj.vel.block(i) ;
			}
		}

		forms.R_visc += W_LUMPING * forms.R ;

		const typename FormMat<WD,WD>::SymType IP = activeProj.vel.Identity() - activeProj.vel ;
		forms.R_visc = activeProj.vel * forms.R_visc * activeProj.vel  + IP ;

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

		// tau:D(w)
		forms.G.clear();
		forms.G.setRows( n );
		forms.G.setCols( m );
		forms.G.cloneIndex( builder.index() ) ;
		forms.G.setBlocksToZero() ;

		// tau:W(w)
		forms.J.clear();
		forms.J.setRows( n );
		forms.J.setCols( m );
		forms.J.cloneIndex( builder.index() ) ;
		forms.J.setBlocksToZero() ;

		builder.integrate_particle( particles, [&]( Index, Scalar w, const D_Itp& l_itp, const D_Dcdx& , const P_Itp& r_itp, const P_Dcdx& r_dc_dx )
		{
			Builder::addTauDu( forms.G, w, l_itp, r_itp, r_dc_dx, dualNodes.indices, fullIndices ) ;
			Builder::addTauWu( forms.J, w, l_itp, r_itp, r_dc_dx, dualNodes.indices, fullIndices ) ;
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

		// tau:D(u)
		forms.H.clear();
		forms.H.setRows( n );
		forms.H.setCols( ma );
		forms.H.cloneIndex( builder.index() ) ;
		forms.H.setBlocksToZero() ;

		builder.integrate_particle( particles, [&]( Index, Scalar w, const D_Itp& l_itp, const D_Dcdx& , const P_Itp& r_itp, const P_Dcdx& r_dc_dx )
		{
			Builder::addTauDu( forms.H, w*sStk, l_itp, r_itp, r_dc_dx, dualNodes.indices, primalNodes.indices ) ;
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

	particles.integratePrimal( activeCells, intPhiPrimal, intPhiVel ) ;
	particles.integrateDual( intPhiDual, intPhiInertia, intPhiOrient, intPhiCohesion ) ;

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
