#include "PhaseSolver.hh"

#include "Phase.hh"
#include "DynParticles.hh"

#include "FormBuilder.hh"
#include "FormBuilder.impl.hh"

#include "utils/Log.hh"

#include <bogus/Core/Block.impl.hpp>

namespace d6 {

const Index PhaseSolver::Active::s_Inactive = -1 ;

struct PhaseMatrices
{
	DynVec M_lumped ;
	typename FormMat<3,3>::SymType A ;

	typename FormMat<3,3>::Type B ;

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
	mats.Pvel.setCols( m );
	mats.Pvel.reserve( m );

	mats.Pstress.setRows( m );
	mats.Pstress.setCols( m );
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

void PhaseSolver::assembleMatrices(const Config &config, const ScalarField &phiInt,
								   PhaseMatrices& mats ) const
{
	typedef const typename MeshType::Location& Loc ;
	typedef const typename MeshType::Interpolation& Itp ;
	typedef const typename MeshType::Derivatives& Dcdx ;

	const Index m = m_phaseNodes.count() ;

	computeProjectors( mats ) ;

	// M_lumped
	mats.M_lumped.resize( 3*m ) ;
	DynMat3::MapType M_map( mats.M_lumped.data(), 3, m ) ;
	M_map.row( 0 ) = phiInt.flatten() ;
	M_map.row( 1 ) = phiInt.flatten() ;
	M_map.row( 2 ) = phiInt.flatten() ;

	FormBuilder builder( phiInt.mesh() ) ;
	builder.reset( m );
	builder.addToIndex( m_phaseNodes.cells, m_phaseNodes.indices, m_phaseNodes.indices );
	builder.makeCompressed();

	// A
	//FIXME
//	mats.A.cloneStructure( builder.index() ) ;
	builder.integrate_qp( m_phaseNodes.cells, [&]( Scalar w, Loc, Itp itp, Dcdx dc_dx )
		{
//			FormBuilder::addDuDv( mats.A, w, itp, dc_dx, m_phaseNodes.indices, m_phaseNodes.indices ) ;
		}
	);

}

void PhaseSolver::step( const Config &config, Phase &phase)
{
	const MeshType& mesh = phase.fraction.mesh() ;

	BoundaryMapper mapper ;
	mesh.make_bc( mapper, m_surfaceNodes ) ;

	// Splat
	ScalarField intPhi    ( mesh ) ;
	VectorField intPhiVel ( mesh ) ;
	ScalarField intPhiInertia( mesh ) ;
	TensorField intPhiOrient ( mesh ) ;
	std::vector< bool > activeCells ;
	m_particles.read( activeCells, intPhi, intPhiVel, intPhiInertia, intPhiOrient ) ;

	// Active nodes
	computeActiveNodes( mesh, activeCells ) ;
	Log::Verbose() << "Active nodes: " << m_phaseNodes.count() << " / " << mesh.nNodes() << std::endl;

	// Matrices
	PhaseMatrices matrices ;
	assembleMatrices( config, intPhi, matrices );

	phase.velocity.set_constant( Vec(.5,0.,-1) );

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

} //d6
