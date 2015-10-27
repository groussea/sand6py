#include "PhaseSolver.hh"

#include "Phase.hh"
#include "DynParticles.hh"

#include "utils/Log.hh"

namespace d6 {

const Index PhaseSolver::Active::s_Inactive = -1 ;

PhaseSolver::PhaseSolver(const DynParticles &particles)
	: m_particles(particles)
{

}


void PhaseSolver::step( const Config &config, Phase &phase)
{
	const MeshType& mesh = phase.fraction.mesh() ;

	BoundaryMapper mapper ;
	mesh.make_bc( mapper, m_surfaceNodes ) ;

	ScalarField intPhi    ( mesh ) ;
	VectorField intPhiVel ( mesh ) ;
	ScalarField intPhiInertia( mesh ) ;
	TensorField intPhiOrient ( mesh ) ;

	std::vector< bool > activeCells ;
	m_particles.read( activeCells, intPhi, intPhiVel, intPhiInertia, intPhiOrient ) ;
	computeActiveNodes( mesh, activeCells ) ;

	Log::Verbose() << "Active nodes: " << m_phaseNodes.count() << " / " << mesh.nNodes() << std::endl;

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
