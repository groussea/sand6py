#include "PhaseSolver.hh"

#include "Phase.hh"
#include "DynParticles.hh"

#include "RigidBody.hh"
#include "RigidBodyData.hh"

#include "FormBuilder.hh"
#include "FormBuilder.impl.hh"

#include "LinearSolver.hh"
#include "solve/Primal.hh"
#include "solve/LCP.hh"

#include "utils/Log.hh"
#include "utils/Config.hh"

#include <bogus/Core/Block.impl.hpp>
#include <bogus/Core/Block.io.hpp>
#include <bogus/Core/Utils/Timer.hpp>

//#define FULL_FEM  // Ignore particles, just solve FEM system

namespace d6 {


struct PhaseMatrices
{
	DynVec S  ;	// Volume aossciated to each node, computed from active cells only

	// Lumped mass matrix, its inverse and a factorization
	typename FormMat<3,3>::SymType M_lumped ;
	typename FormMat<3,3>::SymType M_lumped_inv ;
	typename FormMat<3,3>::SymType M_lumped_inv_sqrt ;

	typename FormMat<3,3>::Type A ; // Mass + Visco ; Could be Symmetric when FormBuilder has sym index

	typename FormMat<6,3>::Type B ; // \phi Tau:D(u)
	typename FormMat<3,3>::Type J ; // \phi Tau:W(u)

	typename FormMat<3,1>::Type C ; // \phi v.grad(p)

	typename FormMat<3,3>::SymType Pvel ;
	typename FormMat<6,6>::SymType Pstress ;

	typename FormMat<6,6>::SymType Aniso ;
};

PhaseSolver::PhaseSolver(const DynParticles &particles)
	: m_particles(particles)
{

}

PhaseSolver::~PhaseSolver()
{
}

void PhaseSolver::computeProjectors(PhaseMatrices& mats,
									const std::vector<RigidBodyData> &rbData) const
{
	const Index m  = m_phaseNodes.count() ;
	const Index mc = nSuppNodes() ;

	mats.Pvel.setRows( m );
	mats.Pvel.reserve( m );

	mats.Pstress.setRows( m+mc );
	mats.Pstress.reserve( m+mc );

	for( Index i = 0 ; i < m_phaseNodes.count() ; ++i  ) {
		const Index idx = m_phaseNodes.revIndices[ i ] ;

		m_surfaceNodes[idx].   velProj( mats.   Pvel.insertBack( i,i ) ) ;
		m_surfaceNodes[idx].stressProj( mats.Pstress.insertBack( i,i ) ) ;
	}

	for( unsigned k = 0 ; k < rbData.size() ; ++k ) {
		for( Index i = 0 ; i < rbData[k].nodes.count() ; ++i  ) {
			const Index   j = rbData[k].nodes.offset + i ;
			const Index idx = rbData[k].nodes.revIndices[ i ] ;

			// Ignore RB-boundary constraints on Dirichlet boundaries
			if( m_surfaceNodes[idx].bc == BoundaryInfo::Stick )
				mats.Pstress.insertBack( j,j ).setZero() ;
			else
				mats.Pstress.insertBack( j,j ).setIdentity() ;
		}
	}

	mats.Pstress.finalize();
	mats.Pvel   .finalize();
}

void PhaseSolver::assembleMatrices(const Config &config, const MeshType &mesh, const ScalarField &phiInt,
								   PhaseMatrices& mats,
								   std::vector< RigidBodyData >& rbData
								   ) const
{
	bogus::Timer timer;

	typedef const typename MeshType::Location& Loc ;
	typedef const typename MeshType::Interpolation& Itp ;
	typedef const typename MeshType::Derivatives& Dcdx ;

	const Index m  = m_phaseNodes.count() ;
	const Index mc = nSuppNodes() ;

	const Scalar dyn_regul = 1.e-8 ;
	const Scalar con_regul = 1.e-8 ;

	computeProjectors( mats, rbData ) ;

	// Lumped mass matrix
	{
		mats.M_lumped.setRows( m );
		mats.M_lumped.setIdentity() ;
		mats.M_lumped_inv.setRows( m );
		mats.M_lumped_inv.setIdentity() ;
		mats.M_lumped_inv_sqrt.setRows( m );
		mats.M_lumped_inv_sqrt.setIdentity() ;

	#pragma omp parallel for
		for( Index i = 0 ; i < m ; ++i ) {
			mats.M_lumped.block( i ) *= phiInt[ m_phaseNodes.revIndices[i] ] ;
		}
	}

	// Other bilinear forms
	{
		timer.reset() ;

		FormBuilder builder( mesh ) ;
		builder.reset( m );
		builder.addToIndex( m_phaseNodes.cells, m_phaseNodes.indices, m_phaseNodes.indices );
		builder.makeCompressed();

		Log::Debug() << "Index computation: " << timer.elapsed() << std::endl ;

		// S -- used to integrate (phi_max - phi)
		mats.S.setZero( m );

		// A
		mats.A.clear();
		mats.A.setRows( m );
		mats.A.setCols( m );
		mats.A.cloneIndex( builder.index() ) ;
		mats.A.setBlocksToZero() ;
		// J
		mats.J.clear();
		mats.J.setRows( m );
		mats.J.setCols( m );
		mats.J.cloneIndex( builder.index() ) ;
		mats.J.setBlocksToZero() ;

		// C
		if( config.enforceMaxFrac ) {
			mats.C.clear();
			mats.C.setRows( builder.rows() );
			mats.C.setCols( m );
			mats.C.cloneIndex( builder.index() ) ;
			mats.C.setBlocksToZero() ;
		}

		builder.addRows(mc) ;

		// B
		mats.B.clear();
		mats.B.setRows( builder.rows() );
		mats.B.setCols( m );
		mats.B.cloneIndex( builder.index() ) ;
		mats.B.setBlocksToZero() ;

		timer.reset() ;
		builder.integrate_qp( m_phaseNodes.cells, [&]( Scalar w, Loc, Itp itp, Dcdx dc_dx )
			{
				FormBuilder:: addDuDv( mats.A, w, itp, dc_dx, m_phaseNodes.indices, m_phaseNodes.indices ) ;
				for( Index k = 0 ; k < MeshType::NQ ; ++k ) {
					mats.S[ m_phaseNodes.indices[itp.nodes[k]] ] += w / MeshType::NQ ;
				}
			}
		);
		Log::Debug() << "Integrate grid: " << timer.elapsed() << std::endl ;


		timer.reset() ;
		builder.integrate_particle( m_particles.geo(), [&]( Index, Scalar w, Loc, Itp itp, Dcdx dc_dx )
			{
				FormBuilder::addTauDu( mats.B, w, itp, dc_dx, m_phaseNodes.indices, m_phaseNodes.indices ) ;
				FormBuilder::addTauWu( mats.J, w, itp, dc_dx, m_phaseNodes.indices, m_phaseNodes.indices ) ;
				if( config.enforceMaxFrac ) {
					FormBuilder::addVDp  ( mats.C, w, itp, dc_dx, m_phaseNodes.indices, m_phaseNodes.indices ) ;
				}
			}
		);
		Log::Debug() << "Integrate particle: " << timer.elapsed() << std::endl ;

	}

	// Rigid bodies
	timer.reset() ;
#pragma omp parallel for if( rbData.size() > 1)
	for( unsigned k = 0 ; k < rbData.size() ; ++k )
	{
		rbData[k].assemble_matrices( m_phaseNodes, m+mc ) ;
	}
	Log::Debug() << "Integrate rbs: " << timer.elapsed() << std::endl ;


	// A = mass + viscosity
#ifndef FULL_FEM
	mats.M_lumped *= config.volMass * config.inv_dt() ;

	mats.A *= 2 * config.viscosity ;
	mats.A += mats.M_lumped ;
#else
	(void) config ;
#endif


	// Projections
	const typename FormMat<3,3>::SymType IP = mats.Pvel.Identity() - mats.Pvel ;
	mats.A = mats.Pvel * ( mats.A * mats.Pvel ) + IP ;

#pragma omp parallel for
	for( Index i = 0 ; i < m ; ++i ) {
		const Scalar m = mats.M_lumped.block(i).trace() / 3 ;
		mats.M_lumped         .block(i) = mats.Pvel.block(i) * m
				+ Mat::Identity() - mats.Pvel.block(i) ;
		mats.M_lumped_inv     .block(i) = mats.Pvel.block(i) * 1./(m + dyn_regul )
				+ Mat::Identity() - mats.Pvel.block(i) ;
		mats.M_lumped_inv_sqrt.block(i) = mats.Pvel.block(i) * 1./std::sqrt( m + con_regul ) ;
	}


}

void PhaseSolver::computeAnisotropy(const DynVec &orientation, const Config& config,
									PhaseMatrices &matrices ) const
{

	matrices.Aniso = matrices.Pstress.Identity() ;

	if( config.anisotropy <= 0 )
		return ;

	const Index m  = m_phaseNodes.count() ;

#pragma omp parallel for
	for( Index i = 0 ; i < m ; ++i ) {

		Mat ori ;
		tensor_view( Segmenter<6>::segment( orientation, i ) ).get( ori ) ;

		ori = (1. - config.anisotropy) * Mat::Identity() + 3 * config.anisotropy * ori ;

		compute_anisotropy_matrix( ori, matrices.Aniso.block(i) );

	}

}

void PhaseSolver::step(const Config &config, Phase &phase,
					   std::vector< RigidBody   >& rigidBodies,
					   std::vector<TensorField > &rbStresses)
{
	bogus::Timer timer ;

	const MeshType& mesh = phase.fraction.mesh() ;

	PhaseMatrices matrices ;
	std::vector< RigidBodyData > rbData ;

	ScalarField intPhi ( mesh ) ;
	DynVec phiu_int ;
	DynVec cohesion ;
	DynVec inertia  ;

	{

		m_surfaceNodes.clear();
		m_surfaceNodes.resize( mesh.nNodes() );

		// Compute volumes of cells
		ScalarField volumes(mesh) ; volumes.set_zero();
		for( typename MeshType::CellIterator it = mesh.cellBegin() ; it != mesh.cellEnd() ; ++it ) {
			typename MeshType::CellGeo geo ;
			typename MeshType::NodeList nodes ;

			mesh.get_geo( *it, geo );
			const Scalar vol = geo.volume() / MeshType::NV ;

			mesh.list_nodes( *it, nodes );
			for( Index k = 0 ; k < MeshType::NV ; ++k ) {
				volumes[nodes[k]] += vol ;
			}
		}

		// Rasterize particles
		VectorField intPhiVel     ( mesh ) ;
		ScalarField intPhiInertia ( mesh ) ;
		ScalarField intPhiCohesion( mesh ) ;
		TensorField intPhiOrient  ( mesh ) ;
		std::vector< bool > activeCells ;


	#if defined(FULL_FEM)
		intPhi.set_constant( 1. ) ;
		intPhiVel.set_zero() ;
		intPhiInertia.set_zero() ;
		intPhiCohesion.set_zero() ;
		intPhiOrient.set_zero() ;
		activeCells.assign( activeCells.size(), true ) ;
	#else
		m_particles.read( activeCells, intPhi, intPhiVel, intPhiInertia, intPhiOrient, intPhiCohesion ) ;
	#endif

		// Compute phi and grad_phi
		phase.fraction.flatten() = intPhi.flatten() ;
		phase.fraction.divide_by( volumes ) ;
		computeGradPhi( phase.fraction, volumes, phase.grad_phi ) ;

		// Active nodes
		computeActiveNodes( activeCells, phase.grad_phi ) ;
		Log::Verbose() << "Active nodes: " << m_phaseNodes.count() << " / " << mesh.nNodes() << std::endl;

		//Rigid bodies and frictional boundaries
		computeActiveBodies( rigidBodies, rbStresses, rbData );

		// Matrices
		mesh.make_bc( StrBoundaryMapper( config.boundary ), m_surfaceNodes ) ;
		assembleMatrices( config, mesh, intPhi, matrices, rbData );

		Log::Debug() << "Matrices assembled  at " << timer.elapsed() << std::endl ;


		m_phaseNodes.field2var( intPhiVel, phiu_int ) ;

		// Cohesion, inertia, orientation
		DynVec orientation ;
		intPhiCohesion.divide_by_positive( intPhi ) ;
		intPhiInertia .divide_by_positive( intPhi ) ;
		intPhiOrient  .divide_by_positive( intPhi ) ;

		m_phaseNodes.field2var( intPhiCohesion, cohesion ) ;
		m_phaseNodes.field2var( intPhiInertia , inertia  ) ;
		m_phaseNodes.field2var( intPhiOrient  , orientation  ) ;

		computeAnisotropy( orientation, config, matrices );
	}

	{

		// Compute fraction of grains
		DynVec fraction ;
		m_phaseNodes.field2var( phase.fraction, fraction ) ;


		// Compute rhs of momentum conservation -- gravity + u(t)
		DynVec rhs ;
		{
			DynVec forces ;

			//Gravity
			{
				VectorField gravity ( mesh ) ;
				gravity.set_constant( config.gravity );
				DynVec grav ; m_phaseNodes.field2var( gravity, grav ) ;
				forces = config.dt() * matrices.M_lumped * grav ;
			}

			// Inertia
#ifndef  FULL_FEM
			forces += config.inv_dt() * config.volMass * phiu_int ;
#endif

			rhs = matrices.Pvel * forces ;
		}

		// Solve unconstrained momentum equation
		DynVec u = matrices.M_lumped_inv * rhs ;
		solveSDP( matrices.A, matrices.M_lumped_inv, rhs, u ) ;

		Log::Debug() << "Linear solve at " << timer.elapsed() << std::endl ;

		// Maximum fraction projection
		if( config.enforceMaxFrac){
			Eigen::VectorXd depl ;
			enforceMaxFrac( config, matrices, rbData, fraction, depl );

			m_phaseNodes.var2field( depl, phase.geo_proj ) ; //Purely geometric
			// u += depl/config.dt() ; // Includes proj into inertia
		}

		// Friction solve
		{

			solveComplementarity( config, matrices, rbData, fraction, cohesion, inertia, u, phase );
			Log::Debug() << "Complementarity solve at " << timer.elapsed() << std::endl ;
		}

		m_phaseNodes.var2field( u, phase.velocity ) ;

		{
			// Velocities gradient
			DynVec int_phiDu = .5 * matrices.Pstress * DynVec( matrices.B * u ) ;
			m_phaseNodes.var2field( int_phiDu, phase.sym_grad ) ;
			phase.sym_grad.divide_by_positive( intPhi ) ;

			DynVec int_phiWu = .5 * matrices.J * u ;
			m_phaseNodes.var2field( int_phiWu, phase.spi_grad ) ;
			phase.spi_grad.divide_by_positive( intPhi ) ;
		}
	}

	Log::Debug() << "Max vel: " << phase.velocity.max_abs() << std::endl ;
}

void PhaseSolver::computeGradPhi( const ScalarField& fraction, const ScalarField& volumes, VectorField &grad_phi ) const
{
	const MeshType &mesh = grad_phi.mesh() ;

	grad_phi.set_zero() ;

	FormBuilder builder(mesh) ;

	typedef const typename MeshType::Location& Loc ;
	typedef const typename MeshType::Interpolation& Itp ;
	typedef const typename MeshType::Derivatives& Dcdx ;

	builder.integrate_qp( mesh.cellBegin(), mesh.cellEnd(),[&]( Scalar w, Loc , Itp itp, Dcdx dc_dx )
	{
		for( Index j = 0 ; j < MeshType::NV ; ++j ) {
			for( Index k = 0 ; k < MeshType::NV ; ++k ) {
				grad_phi[ itp.nodes[j] ] += w * dc_dx.row(k) * itp.coeffs[k] * fraction[ itp.nodes[k] ] ;
			}
		}
	} ) ;

	grad_phi.divide_by( volumes ) ;
}

void PhaseSolver::computeActiveNodes(const std::vector<bool> &activeCells ,
									 const VectorField &grad_phi )
{
	const MeshType &mesh = grad_phi.mesh() ;

	m_phaseNodes.reset( mesh.nNodes() );

	std::vector< int > activeNodes( mesh.nNodes(), 0 ) ;

	Eigen::Matrix< Scalar, 3, Eigen::Dynamic > vecs( 3, mesh.nNodes() ) ;
	vecs.setZero() ;

	for( typename MeshType::CellIterator it = mesh.cellBegin() ; it != mesh.cellEnd() ; ++it )
	{
		if (!activeCells[ it.index() ] ) continue ;

		m_phaseNodes.cells.push_back( *it );

		typename MeshType::NodeList nodes ;
		mesh.list_nodes( *it, nodes );

		for( int k = 0 ; k < MeshType::NV ; ++ k ) {
			++activeNodes[ nodes[k] ] ;
		}


	}

	for( size_t i = 0 ; i < activeNodes.size() ; ++i ) {
		if( activeNodes[i] > 0 ) {

			m_phaseNodes.indices[i] = m_phaseNodes.nNodes++ ;

			m_surfaceNodes[i].bc = BoundaryInfo::Interior ;

			if( activeNodes[i] < mesh.nAdjacent(i) &&
				grad_phi[i].squaredNorm() > 1.e-16) {

				m_surfaceNodes[i].bc = BoundaryInfo::Free ;
				m_surfaceNodes[i].normal = - grad_phi[ i ].normalized() ;

			}
		}
	}

	m_phaseNodes.computeRevIndices();
}

void PhaseSolver::computeActiveBodies( std::vector<RigidBody> &rigidBodies,
							  std::vector<TensorField> &rbStresses,
							  std::vector< RigidBodyData >& rbData )
{
	rbData.clear();
	for( unsigned i = 0 ; i < rigidBodies.size() ; ++i ) {
		rbData.emplace_back( rigidBodies[i], rbStresses[i] );
	}

#pragma omp parallel for
	for( unsigned i = 0 ; i < rigidBodies.size() ; ++i ) {
		rbData[i].compute_active( m_phaseNodes, m_surfaceNodes ) ;
		rbData[i].nodes.computeRevIndices() ;
	}

	m_totRbNodes = 0 ;
	for( unsigned i = 0 ; i < rigidBodies.size() ; ++i )
	{
		rbData[i].nodes.setOffset( m_totRbNodes + m_phaseNodes.count() ) ;
		m_totRbNodes += rbData[i].nodes.count() ;
	}
	Log::Debug() << "Tot coupling nodes: " << nSuppNodes() << std::endl ;



}

void PhaseSolver::solveComplementarity(const Config &c, const PhaseMatrices &matrices,
									   std::vector<RigidBodyData> &rbData,
									   const DynVec &fraction, const DynVec &cohesion, const DynVec &inertia,
									   DynVec &u, Phase& phase ) const
{
	PrimalData	data ;
	data.H = matrices.Aniso * ( matrices.Pstress * ( matrices.B * matrices.M_lumped_inv_sqrt ) ) ;

	data.mu.resize( data.n() ) ;

	// Inertia, mu(I) = \delta_mu * (1./ (1 + I0/I) ), I = dp * sqrt( rho ) * inertia, inertia = |D(U)|/sqrt(p)
	const Scalar I0bar = c.I0 / ( c.grainDiameter * std::sqrt( c.volMass )) ;
	data.mu.segment(0,fraction.rows()).array() = c.mu + c.delta_mu / ( 1. + I0bar / inertia.array().max(1.e-12) ) ;

	//Cohesion : add \grad{ c phi } to rhs
	{
		const Scalar cohe_start = 0.999 * c.phiMax ;
		const DynArr contact_zone = ( ( fraction.array() - cohe_start )/ (c.phiMax - cohe_start) ).max(0).min(1) ;

		DynVec cohe_s  = DynVec::Zero( data.H.rows() ) ;
		component< 6 >( cohe_s, 0 ).head(cohesion.rows()).array() =
				c.cohesion * cohesion.array() * contact_zone ;

		u -= matrices.M_lumped_inv_sqrt * DynVec( data.H.transpose() * cohe_s ) ;
	}

	data.w = matrices.Aniso * DynVec( matrices.Pstress * DynVec( matrices.B * u ) ) ;

	data.jacobians.reserve( rbData.size() ) ;
	data.inv_inertia_matrices.resize( 6, 6*rbData.size() ) ;
	std::vector< unsigned > coupledRbIndices ;

	DynVec totFraction = fraction ;

	// Handle rigid bodies jacobians
	for( unsigned k = 0 ; k < rbData.size() ; ++k ) {
		RigidBodyData& rb = rbData[k] ;

		if( rb.nodes.count() == 0 )
			continue ;

		typename FormMat<6,3>::Type J =	matrices.Aniso * ( matrices.Pstress * ( rb.jacobian ) ) ;

		data.H -= J * matrices.M_lumped_inv_sqrt ;

		const DynVec delta_u = u - rb.projection.transpose() * rb.rb.velocities() ;

		data.w -= J * delta_u  ;

		for( Index i = 0 ; i < rb.nodes.count() ; ++i ) {
			totFraction( m_phaseNodes.indices[ rb.nodes.revIndices[i] ] ) += rb.fraction[i] ;
		}

		data.mu.setConstant( rb.nodes.offset, rb.nodes.count(), c.muRigid ) ;

		Mat66 inv_inertia ;
		rb.rb.inv_inertia( inv_inertia ) ;
		if( inv_inertia.squaredNorm() < 1.e-16 )
			continue ;

		coupledRbIndices.push_back( k ) ;
		data.inv_inertia_matrices.block<6,6>( 0, 6*data.jacobians.size() ) = inv_inertia * c.dt() ;
		data.jacobians.emplace_back( J * rb.projection.transpose() );
	}


	// Compressability
	{
		const DynVec q = ( ( c.phiMax - totFraction.array() )
						   * s_sqrt_23 * matrices.S.array()    // 1/d * Tr \tau = \sqrt{2}{d} \tau_N
						   * c.inv_dt()
						   ).max( 0 ) ;

		component< 6 >( data.w, 0 ).head(q.rows()) += q ;
	}

	DynVec x( data.w.rows() ), y( data.w.rows() ) ;

	// Warm-start stresses
	m_phaseNodes.field2var( phase.stresses, x, false ) ;
	for( unsigned k = 0 ; k < rbData.size() ; ++k ) {
		RigidBodyData& rb = rbData[k] ;
		rb.nodes.field2var( rb.stresses, x, false ) ;
	}

	// Proper solving
	Primal( data ).solve( x, y ) ;

	// Output
	u += matrices.M_lumped_inv_sqrt * DynVec( data.H.transpose() * x ) ;

	// Contact forces -- useless, debug only
	{
		const DynVec fcontact = matrices.Pvel * DynVec( matrices.B.transpose()
														* DynVec( matrices.Pstress
																  * DynVec( matrices.Aniso * x ) ) );
		m_phaseNodes.var2field( fcontact, phase.fcontact ) ;
	}

	// Save stresses for warm-starting next step
	m_phaseNodes.var2field( x, phase.stresses ) ;
	for( unsigned k = 0 ; k < rbData.size() ; ++k ) {
		RigidBodyData& rb = rbData[k] ;
		rb.nodes.var2field( x, rb.stresses ) ;
	}

	// Add contact forces to rigid bodies
	for( unsigned k = 0 ; k < coupledRbIndices.size() ; ++k ) {
		RigidBodyData& rb = rbData[ coupledRbIndices[k] ] ;
		const Vec6 forces = data.jacobians[k].transpose() * x ;
		rb.rb.integrate_forces( c.dt(), forces );
	}

}

void PhaseSolver::enforceMaxFrac(const Config &c, const PhaseMatrices &matrices,
									   std::vector<RigidBodyData> &rbData,
									   const DynVec &fraction,
									   DynVec& depl ) const
{

	LCPData data ;

	data.H = matrices.Pvel * matrices.C ;
	data.w.setZero( data.n() );

	DynVec totFraction = fraction ;

	for( unsigned k = 0 ; k < rbData.size() ; ++k ) {
		RigidBodyData& rb = rbData[k] ;

		if( rb.nodes.count() == 0 )
			continue ;

		for( Index i = 0 ; i < rb.nodes.count() ; ++i ) {
			Scalar& frac = totFraction( m_phaseNodes.indices[ rb.nodes.revIndices[i] ] ) ;
			frac = std::min( std::max( 1., frac ), frac + rb.fraction[i] ) ;
		}
	}

	data.w.segment(0, totFraction.rows()) = ( c.phiMax - totFraction.array() )
			* matrices.S.array()    // 1/d * Tr \tau = \sqrt{2}{d} \tau_N
			;

	// Apply pressure projection
#pragma omp parallel for
	for( Index i = 0 ; i < data.H.rowsOfBlocks() ; ++i ) {
		for( LCPData::HType::InnerIterator it = data.H.innerIterator(i) ; it ; ++it ) {
			data.H.block( it.ptr() ) *= matrices.Pstress.block( it.inner() )(0,0) ;
		}
	}
	for( Index i = 0 ; i < data.n() ; ++i ) {
		data.w(i) *= matrices.Pstress.block( i )(0,0) ;
	}

	DynVec x = DynVec::Zero( data.n() ) ;
	LCP lcp( data ) ;
	lcp.solve( x ) ;

	depl = - data.H * x ;

}


} //d6
