#include "PhaseSolver.hh"

#include "Phase.hh"
#include "DynParticles.hh"
#include "Config.hh"

#include "FormBuilder.hh"
#include "FormBuilder.impl.hh"

#include "LinearSolver.hh"

#include "solve/Primal.hh"

#include "utils/Log.hh"

#include <bogus/Core/Block.impl.hpp>
#include <bogus/Core/Block.io.hpp>
#include <bogus/Core/Utils/Timer.hpp>

//#define FULL_FEM

namespace d6 {

const Index PhaseSolver::Active::s_Inactive = -1 ;

template < typename Derived >
void PhaseSolver::Active::field2var( const FieldBase<Derived> &field, DynVec & var ) const
{
	const Index n = field.size() ;
	constexpr Index D = FieldBase<Derived>::D ;

	var.resize( count() * D ) ;
#pragma omp parallel for
	for( Index i = 0 ; i < n ; ++ i) {
		const Index idx = indices[ i ] ;
		if( idx != Active::s_Inactive ) {
			Segmenter<D>::segment( var, idx ) = field[ i ] ;
		}
	}
}

template < typename Derived >
void PhaseSolver::Active::var2field( const DynVec & var,  FieldBase<Derived> &field ) const
{
	const Index n = field.mesh().nNodes() ;
	constexpr Index D = FieldBase<Derived>::D ;

	field.set_zero();

#pragma omp parallel for
	for( Index i = 0 ; i < n ; ++ i) {
		const Index idx = indices[ i ] ;
		if( idx != Active::s_Inactive ) {
			field[ i ] = Segmenter<D>::segment( var, idx ) ;
		}
	}
}

struct PhaseMatrices
{
	DynVec S  ;
	DynVec S6 ;

	typename FormMat<3,3>::SymType M_lumped ;
	typename FormMat<3,3>::SymType M_lumped_inv ;
	typename FormMat<3,3>::SymType M_lumped_inv_sqrt ;

	typename FormMat<3,3>::Type A ; //Could be Symmetric when FormBuilder has sym index

	typename FormMat<6,3>::Type B ;
	typename FormMat<3,3>::Type J ;

	typename FormMat<3,3>::SymType Pvel ;
	typename FormMat<6,6>::SymType Pstress ;

};

PhaseSolver::PhaseSolver(const DynParticles &particles)
	: m_particles(particles)
{

}

void PhaseSolver::computeProjectors( PhaseMatrices& mats ) const
{
	const Index m = m_phaseNodes.count() ;

	mats.Pvel.setRows( m );
	mats.Pvel.reserve( m );

	mats.Pstress.setRows( m );
	mats.Pstress.reserve( m );

	for( size_t i = 0 ; i < m_surfaceNodes.size() ; ++i  ) {
		const Index idx = m_phaseNodes.indices[ i ] ;
		if( idx != Active::s_Inactive ) {
			m_surfaceNodes[i].   velProj( mats.   Pvel.insertBack( idx,idx ) ) ;
			m_surfaceNodes[i].stressProj( mats.Pstress.insertBack( idx,idx ) ) ;
		}
	}

	mats.Pstress.finalize();
	mats.Pvel   .finalize();
}

void PhaseSolver::assembleMatrices(const Config &config, const MeshType &mesh, const DynVec &phiInt,
								   PhaseMatrices& mats ) const
{
	typedef const typename MeshType::Location& Loc ;
	typedef const typename MeshType::Interpolation& Itp ;
	typedef const typename MeshType::Derivatives& Dcdx ;

	const Index m = m_phaseNodes.count() ;

	const Scalar regul = 1.e-8 ;

	computeProjectors( mats ) ;

	// M_lumped
//	mats.M_lumped.resize( 3*m ) ;
//	DynMat3::MapType M_map( mats.M_lumped.data(), 3, m ) ;
//	M_map.row( 0 ) = phiInt.flatten() ;
//	M_map.row( 1 ) = phiInt.flatten() ;
//	M_map.row( 2 ) = phiInt.flatten() ;

	mats.M_lumped.setRows( m );
	mats.M_lumped.setIdentity() ;
	mats.M_lumped_inv.setRows( m );
	mats.M_lumped_inv.setIdentity() ;
	mats.M_lumped_inv_sqrt.setRows( m );
	mats.M_lumped_inv_sqrt.setIdentity() ;

#pragma omp parallel for
	for( Index i = 0 ; i < m ; ++i ) {
		mats.M_lumped.block( i ) *= phiInt[i] ;
	}

	FormBuilder builder( mesh ) ;
	builder.reset( m );
	builder.addToIndex( m_phaseNodes.cells, m_phaseNodes.indices, m_phaseNodes.indices );
	builder.makeCompressed();

	// S
	mats.S .resize(   m ); mats.S .setZero();
	mats.S6.resize( 6*m ); mats.S6.setZero();

	// A
	mats.A.clear();
	mats.A.setRows( m );
	mats.A.setCols( m );
	mats.A.cloneIndex( builder.index() ) ;
	mats.A.setBlocksToZero() ;

	// B
	mats.B.clear();
	mats.B.setRows( m );
	mats.B.setCols( m );
	mats.B.cloneIndex( builder.index() ) ;
	mats.B.setBlocksToZero() ;

	// J
	mats.J.clear();
	mats.J.setRows( m );
	mats.J.setCols( m );
	mats.J.cloneIndex( builder.index() ) ;
	mats.J.setBlocksToZero() ;

	builder.integrate_qp( m_phaseNodes.cells, [&]( Scalar w, Loc, Itp itp, Dcdx dc_dx )
		{
			FormBuilder:: addDuDv( mats.A, w, itp, dc_dx, m_phaseNodes.indices, m_phaseNodes.indices ) ;
			for( Index k = 0 ; k < MeshType::NQ ; ++k ) {
										mats.S[ m_phaseNodes.indices[itp.nodes[k]] ] += w / MeshType::NQ ;
				Segmenter<6>::segment( mats.S6, m_phaseNodes.indices[itp.nodes[k]] ) += Vec6::Constant( w / MeshType::NQ ) ;
			}

		}
	);

	builder.integrate_particle( m_particles.geo(), [&]( Scalar w, Loc, Itp itp, Dcdx dc_dx )
		{
			FormBuilder::addTauDu( mats.B, w, itp, dc_dx, m_phaseNodes.indices, m_phaseNodes.indices ) ;
			FormBuilder::addTauWu( mats.J, w, itp, dc_dx, m_phaseNodes.indices, m_phaseNodes.indices ) ;
		}
	);

	/////////////////
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

//#pragma omp parallel for
	for( Index i = 0 ; i < m ; ++i ) {
		const Scalar m = mats.M_lumped.block(i).trace() / 3 ;
		mats.M_lumped         .block(i) = mats.Pvel.block(i) * m
				+ Mat::Identity() - mats.Pvel.block(i) ;
		mats.M_lumped_inv     .block(i) = mats.Pvel.block(i) * 1./( regul + m )
				+ Mat::Identity() - mats.Pvel.block(i) ;
		mats.M_lumped_inv_sqrt.block(i) = mats.Pvel.block(i) * 1./std::sqrt( regul + m )
				+ Mat::Identity() - mats.Pvel.block(i) ;
	}

}

struct CuveBoundary: public BoundaryMapper
{
	virtual BoundaryInfo::Bc operator() ( const std::string & domain ) const
	{
		if( domain == "top")
			return BoundaryInfo::Normal ;
		if( domain == "bottom")
			return BoundaryInfo::Stick ;

		return BoundaryInfo::Slip ;
	}

};

void PhaseSolver::step( const Config &config, Phase &phase)
{
	bogus::Timer timer ;

	const MeshType& mesh = phase.fraction.mesh() ;
	m_surfaceNodes.resize( mesh.nNodes() );

	VectorField gravity ( mesh ) ;
	gravity.set_constant( config.gravity );

	// Splat
	ScalarField intPhi    ( mesh ) ;
	VectorField intPhiVel ( mesh ) ;
	ScalarField intPhiInertia( mesh ) ;
	TensorField intPhiOrient ( mesh ) ;
	std::vector< bool > activeCells ;
	m_particles.read( activeCells, intPhi, intPhiVel, intPhiInertia, intPhiOrient ) ;

#ifdef FULL_FEM
	activeCells.assign( activeCells.size(), true ) ;
#endif

	// Active nodes
	computeActiveNodes( mesh, activeCells ) ;
	Log::Verbose() << "Active nodes: " << m_phaseNodes.count() << " / " << mesh.nNodes() << std::endl;

	DynVec phi_int  ; m_phaseNodes.field2var( intPhi, phi_int ) ;
	DynVec phiu_int ; m_phaseNodes.field2var( intPhiVel, phiu_int ) ;

#ifdef FULL_FEM
	phi_int.setOnes() ;
	phiu_int.setZero() ;
#endif

	// Matrices

	BoundaryMapper mapper ;
//	CuveBoundary mapper ;
	mesh.make_bc( mapper, m_surfaceNodes ) ;
	PhaseMatrices matrices ;
	assembleMatrices( config, mesh, phi_int, matrices );

	Log::Debug() << "Matrices assembled  at " << timer.elapsed() << std::endl ;
	{
#ifdef FULL_FEM
		phi_int = matrices.S ;
#endif
		DynVec fraction = phi_int.array() / matrices.S.array() ;
		m_phaseNodes.var2field( fraction, phase.fraction ) ;

		DynVec rhs ;
		{
			DynVec grav ;	m_phaseNodes.field2var( gravity, grav ) ;
			rhs = matrices.Pvel * DynVec( config.inv_dt() * config.volMass * phiu_int
										  + config.dt() * matrices.M_lumped * grav ) ;
		}

		DynVec u = matrices.M_lumped_inv * rhs ;
		solveSDP( matrices.A, matrices.M_lumped_inv, rhs, u ) ;

		Log::Debug() << "Linear solve at " << timer.elapsed() << std::endl ;

		solveComplementarity( config, matrices, fraction, u, phase );

		Log::Debug() << "Complementarity solve at " << timer.elapsed() << std::endl ;

		m_phaseNodes.var2field( u, phase.velocity ) ;

		{
			// Velocities gradient
			DynVec int_phiDu = .5 * matrices.Pstress * DynVec( matrices.B * u ) ;
			m_phaseNodes.var2field( int_phiDu, phase.sym_grad ) ;
			phase.sym_grad.divide_by( intPhi ) ;

			DynVec int_phiWu = .5 * matrices.J * u ;
			m_phaseNodes.var2field( int_phiWu, phase.spi_grad ) ;
			phase.spi_grad.divide_by( intPhi ) ;
		}
	}

	Log::Debug() << "Max vel: " << phase.velocity.max_abs() << std::endl ;
}

void PhaseSolver::computeActiveNodes(
		const MeshType& mesh, const std::vector<bool> &activeCells )
{
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

		typename MeshType::CellGeo cellGeo ;
		mesh.get_geo( *it, cellGeo );

		for( int k = 0 ; k < MeshType::NV ; ++ k ) {
			++activeNodes[ nodes[k] ] ;
			vecs.col( nodes[k] ) += ( cellGeo.vertex(k) - cellGeo.center() ).normalized() ;
		}
	}

	for( size_t i = 0 ; i < activeNodes.size() ; ++i ) {
		if( activeNodes[i] > 0 ) {

			m_phaseNodes.indices[i] = m_phaseNodes.nNodes ;
			++ m_phaseNodes.nNodes ;

			m_surfaceNodes[i].bc = BoundaryInfo::Interior ;

			if( activeNodes[i] < mesh.nAdjacent(i) ) {

				m_surfaceNodes[i].bc = BoundaryInfo::Free ;
				m_surfaceNodes[i].normal = vecs.col( i ) / activeNodes[i] ;

			}
		}
	}

}

void PhaseSolver::solveComplementarity(const Config &c, const PhaseMatrices &matrices,
									   const DynVec &fraction,
									   DynVec &u, Phase& phase ) const
{
	PrimalData	data ;

	data.H = matrices.Pstress * ( matrices.B * matrices.M_lumped_inv_sqrt ) ;
	data.w = matrices.Pstress * DynVec( matrices.B * u ) ;

	data.mu.setConstant( data.n(), c.mu ) ;

	{
		// Compressability
		DynVec q =  c.phiMax * matrices.S ;
		q.array() -= fraction.array() * matrices.S.array();
		q *= c.inv_dt() ;

		q = q.cwiseMax( 0 ) ;

		DynMat6::Map( data.w.data(), 6, data.n() ).row( 0 )	+= q ;
	}

	DynVec x( data.w.rows() ), y( data.w.rows() ) ;
	m_phaseNodes.field2var( phase.stresses, x ) ;

	Primal primal( data ) ;
	primal.solve( x, y ) ;

	u += matrices.M_lumped_inv_sqrt * DynVec( data.H.transpose() * x ) ;

	m_phaseNodes.var2field( x, phase.stresses ) ;

}


} //d6
