#ifndef D6_PHASE_MESHES_HH
#define D6_PHASE_MESHES_HH

#include "PhaseFields.hh"

namespace d6 {

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
