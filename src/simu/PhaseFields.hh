#ifndef D6_PHASE_FIELDS_HH
#define D6_PHASE_FIELDS_HH

#include "geo/geo.fwd.hh"

#include <memory>

#define D6_UNSTRUCTURED_DUAL

namespace d6 {

typedef MeshImpl PrimalMesh ;

typedef Linear<PrimalMesh> PrimalShape ;
//typedef   P2<PrimalMesh> PrimalShape ;

#ifdef D6_UNSTRUCTURED_DUAL
typedef UnstructuredDOFs      DualMesh ;
typedef UnstructuredShapeFunc DualShape ;
#else
typedef MeshImpl   DualMesh ;

typedef DGLinear<  DualMesh> DualShape ;
//typedef DGConstant<  DualMesh> DualShape ;
//typedef   Linear<  DualMesh> DualShape ;
#endif

typedef AbstractScalarField< PrimalShape > PrimalScalarField ;
typedef AbstractVectorField< PrimalShape > PrimalVectorField ;
typedef AbstractTensorField< PrimalShape > PrimalTensorField ;

typedef AbstractScalarField<   DualShape > DualScalarField ;
typedef AbstractTensorField<   DualShape > DualTensorField ;
typedef AbstractSkewTsField<   DualShape > DualSkewTsField ;

typedef PrimalTensorField RBStresses ;

class  DynParticles ;
struct Phase ;
struct Config ;

template <typename PMeshT, typename DMeshT>
struct AbstractPhaseMeshes
{
	std::unique_ptr<PrimalMesh> m_primal ;
	std::unique_ptr<  DualMesh> m_dual  ;

	const PrimalMesh& primal() const { return *m_primal ; }
	const   DualMesh&   dual() const { return   *m_dual ; }

	AbstractPhaseMeshes( std::unique_ptr<PrimalMesh> pMesh, std::unique_ptr<  DualMesh> dMesh )
		: m_primal( std::move(pMesh)),
		  m_dual  ( std::move(dMesh))
	{}

	void adapt( const DynParticles& particles, std::unique_ptr< Phase >& phase ) ;

	template <typename Ar>
	void serialize( Ar &ar, unsigned int ) {
		ar & (*m_primal) ;
		ar & (*m_dual) ;
	}
};

template <typename PMeshT>
struct AbstractPhaseMeshes<PMeshT, PMeshT> {
	std::unique_ptr<PrimalMesh> m_primal ;

	AbstractPhaseMeshes( std::unique_ptr<PrimalMesh> pMesh, std::unique_ptr<  DualMesh>  )
		: m_primal( std::move( pMesh ) )
	{}

	const PrimalMesh& primal() const { return *m_primal ; }
	const   DualMesh&   dual() const { return *m_primal ; }

	void adapt( const DynParticles& particles, std::unique_ptr< Phase >& phase ) ;

	template <typename Ar>
	void serialize( Ar &ar, unsigned int ) {
		ar & (*m_primal) ;
	}
};

typedef AbstractPhaseMeshes< PrimalMesh, DualMesh > PhaseMeshes ;


} //d6


#endif
