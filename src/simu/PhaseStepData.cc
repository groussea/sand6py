#include "PhaseStepData.hh"

#include "Phase.hh"
#include "DynParticles.hh"
#include "RigidBody.hh"
#include "RigidBodyData.hh"

#include "FormBuilder.impl.hh"

#include "geo/MeshImpl.hh"

#include "utils/Config.hh"
#include "utils/Log.hh"

#include <bogus/Core/Block.impl.hpp>
#include <bogus/Core/Utils/Timer.hpp>

#include <utility>

//#define FULL_FEM  // Ignore particles, just solve FEM system

#define INTEGRATE_PARTICLES_SEQUENTIAL

namespace d6 {

void PhaseStepData::computeProjectors(const bool weakStressBC, const std::vector<RigidBodyData> &rbData,
									  Projectors& mats) const
{
	const Index m  = nNodes() ;
	const Index mc = nSuppNodes() ;

	mats.vel.setRows( m );
	mats.vel.reserve( m );

	mats.stress.setRows( m+mc );
	mats.stress.reserve( m+mc );

	for( Index i = 0 ; i < nodes.count() ; ++i  ) {
		const Index idx = nodes.revIndices[ i ] ;

		m_surfaceNodes[idx].velProj( mats.vel.insertBack( i,i ) ) ;

		if( weakStressBC )
			mats.stress.insertBack( i,i ).setIdentity() ;
		else
			m_surfaceNodes[idx].stressProj( mats.stress.insertBack( i,i ) ) ;
	}

	// Additional nodes for frictional boundaries
	for( unsigned k = 0 ; k < rbData.size() ; ++k ) {
		for( Index i = 0 ; i < rbData[k].nodes.count() ; ++i  ) {
			const Index   j = rbData[k].nodes.offset + i ;
			const Index idx = rbData[k].nodes.revIndices[ i ] ;

			// Ignore RB-boundary constraints on Dirichlet boundaries
			if( m_surfaceNodes[idx].bc == BoundaryInfo::Stick )
				mats.stress.insertBack( j,j ).setZero() ;
			else
				mats.stress.insertBack( j,j ).setIdentity() ;
		}
	}

	mats.stress.finalize();
	mats.vel   .finalize();
}

void PhaseStepData::assembleMatrices(
		const Particles &particles,
		const Config &config, const Scalar dt, const MeshType &mesh, const ScalarField &phiInt,
		std::vector< RigidBodyData >&rbData
		)
{
	// FXIME other approxes

	typedef typename Linear<MeshImpl>::Interpolation Itp ;
	typedef typename Linear<MeshImpl>::Derivatives Dcdx ;

	bogus::Timer timer;

	const Index m  = nNodes() ;
	const Index mc = nSuppNodes() ;

	const Scalar mass_regul = 1.e-8 ;

	computeProjectors( config.weakStressBC, rbData, proj ) ;

	// Lumped mass matrix
	{
		forms.M_lumped.setRows( m );
		forms.M_lumped.setIdentity() ;
		forms.M_lumped_inv.setRows( m );
		forms.M_lumped_inv.setIdentity() ;
		forms.M_lumped_inv_sqrt.setRows( m );
		forms.M_lumped_inv_sqrt.setIdentity() ;

	#pragma omp parallel for
		for( Index i = 0 ; i < m ; ++i ) {
			forms.M_lumped.block( i ) *= phiInt[ nodes.revIndices[i] ] ;
		}
	}

	// Other bilinear forms
	{
		timer.reset() ;

		typedef FormBuilder< Linear<MeshImpl>, Linear<MeshImpl> > Builder ;
		Builder builder( mesh.shaped<Linear>(), mesh.shaped<Linear>() ) ;
		builder.reset( m );
		builder.addToIndex< form::Left >( nodes.cells.begin(), nodes.cells.end(), nodes.indices, nodes.indices );
		builder.makeCompressed();

		Log::Debug() << "Index computation: " << timer.elapsed() << std::endl ;

		// A
		forms.A.clear();
		forms.A.setRows( m );
		forms.A.setCols( m );
		forms.A.cloneIndex( builder.index() ) ;
		forms.A.setBlocksToZero() ;
		// J
		forms.J.clear();
		forms.J.setRows( m );
		forms.J.setCols( m );
		forms.J.cloneIndex( builder.index() ) ;
		forms.J.setBlocksToZero() ;

		// C
		if( config.enforceMaxFrac ) {
			forms.C.clear();
			forms.C.setRows( builder.rows() );
			forms.C.setCols( m );
			forms.C.cloneIndex( builder.index() ) ;
			forms.C.setBlocksToZero() ;
		}

		builder.addRows(mc) ;

		// B
		forms.B.clear();
		forms.B.setRows( builder.rows() );
		forms.B.setCols( m );
		forms.B.cloneIndex( builder.index() ) ;
		forms.B.setBlocksToZero() ;

		timer.reset() ;
#ifdef FULL_FEM
		builder.integrate_qp( nodes.cells, [&]( Scalar w, const Loc&, const Itp& itp, const Dcdx& dc_dx )
			{
				Builder:: addDuDv( forms.A, w, itp, dc_dx, nodes.indices, nodes.indices ) ;
			}
		);
		Log::Debug() << "Integrate grid: " << timer.elapsed() << std::endl ;
#endif

		timer.reset() ;

#ifndef INTEGRATE_PARTICLES_SEQUENTIAL
		{
			// Integrating over particles can be slow and not directly parallelizable
			// To regain some parallelism, we first associate particles to nodes,
			// then compute separately each row of he form matrices

			const size_t n = particles.count() ;

			Eigen::Matrix< Scalar, MeshType::NV, Eigen::Dynamic > coeffs ;
			Eigen::Matrix<  Index, MeshType::NV, Eigen::Dynamic > nodeIds  ;
			coeffs  .resize( MeshType::NV, n) ;
			nodeIds.resize( MeshType::NV, n) ;

			std::vector< std::vector< std::pair< size_t, Index > > > nodeParticles ( m ) ;


			Loc loc ;
			Itp itp ;
			Dcdx dc_dx ;

#pragma omp parallel private( loc, itp )
			for ( size_t i = 0 ; i < n ; ++i )
			{
				mesh.locate( particles.centers().col(i), loc );
				mesh.shaped<Linear>().interpolate( loc, itp );
				nodeIds.col(i) = itp.nodes ;
				coeffs .col(i) = itp.coeffs ;
			}

			for ( size_t i = 0 ; i < n ; ++i ) {
				for( Index k = 0 ; k < MeshType::NV ; ++k ) {
					nodeParticles[ nodes.indices[nodeIds(k,i)] ].push_back( std::make_pair(i,k) ) ;
				}
			}

#pragma omp parallel for private( loc, itp, dc_dx )
			for( Index nidx = 0 ; nidx < m ; ++ nidx )
			{
				for( unsigned i = 0 ; i < nodeParticles[nidx].size() ; ++i ) {
					const size_t pid = nodeParticles[nidx][i].first ;
					const Index k0   = nodeParticles[nidx][i].second ;
					mesh.locate( particles.centers().col(pid), loc );
					mesh.shaped<Linear>().get_derivatives( loc, dc_dx );
					itp.nodes  = nodeIds.col(pid) ;
					itp.coeffs = coeffs.col(pid) ;

					const Scalar w = particles.volumes()[pid] ;
					const Scalar m = itp.coeffs[k0] * w ;
					const Dcdx& const_dcdx = dc_dx ;

					FormBuilder::addTauDu( forms.B, m, nidx, itp, dc_dx, nodes.indices ) ;
					FormBuilder::addTauWu( forms.J, m, nidx, itp, dc_dx, nodes.indices ) ;
					FormBuilder:: addDuDv( forms.A, w, nidx, const_dcdx.row(k0), itp, dc_dx, nodes.indices ) ;
					if( config.enforceMaxFrac ) {
						FormBuilder::addVDp  ( forms.C, m, nidx, itp, dc_dx, nodes.indices ) ;
					}

				}
			}


		}
#else
		builder.integrate_particle( particles, [&]( Index, Scalar w, const Itp& l_itp, const Dcdx& l_dc_dx, const Itp& r_itp, const Dcdx& r_dc_dx )
			{
				Builder::addTauDu( forms.B, w, l_itp, r_itp, r_dc_dx, nodes.indices, nodes.indices ) ;
				Builder::addTauWu( forms.J, w, l_itp, r_itp, r_dc_dx, nodes.indices, nodes.indices ) ;
				Builder:: addDuDv( forms.A, w, l_itp, l_dc_dx, r_itp, r_dc_dx, nodes.indices, nodes.indices ) ;
				if( config.enforceMaxFrac ) {
					Builder::addVDp  ( forms.C, w, l_itp, r_itp, r_dc_dx, nodes.indices, nodes.indices ) ;
				}
			}
		);
#endif
		Log::Debug() << "Integrate particle: " << timer.elapsed() << std::endl ;

	}

	// Rigid bodies
	timer.reset() ;
#pragma omp parallel for if( rbData.size() > 1)
	for( unsigned k = 0 ; k < rbData.size() ; ++k )
	{
		rbData[k].assemble_matrices( nodes, m+mc ) ;
	}
	Log::Debug() << "Integrate rbs: " << timer.elapsed() << std::endl ;


	// A = mass + viscosity
#ifndef FULL_FEM
	forms.M_lumped *= config.volMass / dt  ;

	forms.A *= 2 * config.viscosity ;
	forms.A += forms.M_lumped ;
#endif


	// Projections
	const typename FormMat<WD,WD>::SymType IP = proj.vel.Identity() - proj.vel ;
	forms.A = proj.vel * ( forms.A * proj.vel ) + IP ;

#pragma omp parallel for
	for( Index i = 0 ; i < m ; ++i ) {
		const Scalar mass = forms.M_lumped.block(i).trace() / WD ;
		forms.M_lumped         .block(i) = proj.vel.block(i) * mass
				+ Mat::Identity() - proj.vel.block(i) ;
		forms.M_lumped_inv     .block(i) = proj.vel.block(i) * 1./(mass + mass_regul )
				+ Mat::Identity() - proj.vel.block(i) ;
		forms.M_lumped_inv_sqrt.block(i) = proj.vel.block(i) * 1./std::sqrt( mass + mass_regul ) ;
	}


}

void PhaseStepData::computeAnisotropy(const DynVec &orientation, const Config& config,
									 typename FormMat<SD,SD>::SymType &Aniso ) const
{
	// Compute anisotropy matrix from interpolated orientation distributions

	Aniso = proj.stress.Identity() ;

	if( config.anisotropy <= 0 )
		return ;

	const Index m  = nodes.count() ;

#pragma omp parallel for
	for( Index i = 0 ; i < m ; ++i ) {

		Mat ori ;
		tensor_view( Segmenter<SD>::segment( orientation, i ) ).get( ori ) ;

		ori = (1. - config.anisotropy) * Mat::Identity() + WD * config.anisotropy * ori ;

		compute_anisotropy_matrix( ori, Aniso.block(i) );

	}

}

void PhaseStepData::compute(const DynParticles& particles,
		const Config &config, const Scalar dt, Phase &phase,
		std::vector< RigidBody   >& rigidBodies,
		std::vector<TensorField > &rbStresses , std::vector<RigidBodyData> & rbData)
{
	const ScalarField::ShapeFuncType& shape = phase.fraction.shape() ;
	const MeshType& mesh = shape.derived().mesh() ;

	m_surfaceNodes.clear();
	m_surfaceNodes.resize( mesh.nNodes() );

	// Compute volumes of cells
	ScalarField volumes(shape) ;
	shape.compute_volumes( volumes.flatten() );

	// Transfer particles quantities to grid
	ScalarField intPhi ( shape ) ;
	VectorField intPhiVel     ( shape ) ;
	ScalarField intPhiInertia ( shape ) ;
	ScalarField intPhiCohesion( shape ) ;
	TensorField intPhiOrient  ( shape ) ;
	std::vector< bool > activeCells ;

#if defined(FULL_FEM)
	intPhi.set_constant( 1. ) ;
	intPhiVel.set_zero() ;
	intPhiInertia.set_zero() ;
	intPhiCohesion.set_zero() ;
	intPhiOrient.set_zero() ;
	activeCells.assign( activeCells.size(), true ) ;
#else
	particles.read( activeCells, intPhi, intPhiVel, intPhiInertia, intPhiOrient, intPhiCohesion ) ;
#endif

	// Compute phi and grad_phi (for visualization purposes )
	phase.fraction.flatten() = intPhi.flatten() ;
	phase.fraction.divide_by( volumes ) ;
	computeGradPhi( phase.fraction, volumes, phase.grad_phi ) ;

	// Active nodes
	computeActiveNodes( activeCells, phase.grad_phi ) ;
	Log::Verbose() << "Active nodes: " << nodes.count() << " / " << mesh.nNodes() << std::endl;

	//Rigid bodies and frictional boundaries
	computeActiveBodies( rigidBodies, rbStresses, rbData );

	// Bilinear forms matrices
	mesh.make_bc( StrBoundaryMapper( config.boundary ), m_surfaceNodes ) ;
	nodes.field2var( volumes, forms.volumes ) ;
	assembleMatrices( particles.geo(), config, dt, mesh, intPhi, rbData );

	nodes.field2var( phase.fraction, fraction ) ;
	nodes.field2var( intPhiVel, forms.phiu ) ;

	// Cohesion, inertia, orientation
	DynVec orientation ;
	intPhiCohesion.divide_by_positive( intPhi ) ;
	intPhiInertia .divide_by_positive( intPhi ) ;
	intPhiOrient  .divide_by_positive( intPhi ) ;

	nodes.field2var( intPhiCohesion, cohesion ) ;
	nodes.field2var( intPhiInertia , inertia  ) ;
	nodes.field2var( intPhiOrient  , orientation  ) ;

	computeAnisotropy( orientation, config, Aniso );

	// External forces
	VectorField gravity ( shape ) ;
	gravity.set_constant( config.gravity );
	gravity.multiply_by( intPhi ) ;
	gravity.flatten() *= config.volMass ;
	nodes.field2var( gravity, forms.externalForces ) ;
}


void PhaseStepData::computeGradPhi( const ScalarField& fraction, const ScalarField& volumes, VectorField &grad_phi ) const
{
	const ScalarField::ShapeFuncType& shape = fraction.shape() ;
	const MeshType& mesh = shape.derived().mesh() ;

	grad_phi.set_zero() ;

	// FXIME other approxes
	typedef FormBuilder< Linear<MeshImpl>, Linear<MeshImpl> > Builder ;
	Builder builder( mesh.shaped<Linear>(), mesh.shaped<Linear>() ) ;

	typedef const typename Linear<MeshImpl>::Interpolation& Itp ;
	typedef const typename Linear<MeshImpl>::Derivatives& Dcdx ;

	builder.integrate_qp<form::Left>( mesh.cellBegin(), mesh.cellEnd(),
						  [&]( Scalar w, const Vec&, Itp itp, Dcdx dc_dx, Itp , Dcdx )
	{
		for( Index j = 0 ; j < Linear<MeshImpl>::NI ; ++j ) {
			for( Index k = 0 ; k < Linear<MeshImpl>::NI ; ++k ) {
				grad_phi[ itp.nodes[j] ] += w * dc_dx.row(k) * itp.coeffs[k] * fraction[ itp.nodes[k] ] ;
			}
		}
	} ) ;

	grad_phi.divide_by( volumes ) ;
}

void PhaseStepData::computeActiveNodes(const std::vector<bool> &activeCells ,
									 const VectorField &grad_phi )
{
	const VectorField::ShapeFuncImpl& shape = grad_phi.shape() ;
	const MeshType& mesh = shape.derived().mesh() ;

	nodes.reset( mesh.nNodes() );

	std::vector< int > activeNodes( mesh.nNodes(), 0 ) ;

	Eigen::Matrix< Scalar, WD, Eigen::Dynamic > vecs( WD, mesh.nNodes() ) ;
	vecs.setZero() ;

	for( typename MeshType::CellIterator it = mesh.cellBegin() ; it != mesh.cellEnd() ; ++it )
	{
		if (!activeCells[ it.index() ] ) continue ;

		nodes.cells.push_back( *it );

		typename VectorField::ShapeFuncType::NodeList nodes ;
		shape.list_nodes( *it, nodes );

		for( int k = 0 ; k < Linear<MeshImpl>::NI ; ++ k ) {
			++activeNodes[ nodes[k] ] ;
		}


	}

	for( size_t i = 0 ; i < activeNodes.size() ; ++i ) {
		if( activeNodes[i] > 0 ) {

			nodes.indices[i] = nodes.nNodes++ ;

			m_surfaceNodes[i].bc = BoundaryInfo::Interior ;
		}
	}

	nodes.computeRevIndices();
}

void PhaseStepData::computeActiveBodies( std::vector<RigidBody> &rigidBodies,
										 std::vector<TensorField> &rbStresses,
										 std::vector< RigidBodyData > &rbData
										 )
{
	rbData.clear();
	for( unsigned i = 0 ; i < rigidBodies.size() ; ++i ) {
		rbData.emplace_back( rigidBodies[i], rbStresses[i] );
	}

#pragma omp parallel for
	for( unsigned i = 0 ; i < rigidBodies.size() ; ++i ) {
		rbData[i].compute_active( nodes, m_surfaceNodes ) ;
		rbData[i].nodes.computeRevIndices() ;
	}

	m_totRbNodes = 0 ;
	for( unsigned i = 0 ; i < rigidBodies.size() ; ++i )
	{
		rbData[i].nodes.setOffset( m_totRbNodes + nodes.count() ) ;
		m_totRbNodes += rbData[i].nodes.count() ;
	}
	Log::Debug() << "Tot coupling nodes: " << nSuppNodes() << std::endl ;



}


} //d6
