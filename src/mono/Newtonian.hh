#ifndef D6_NEWTONIAN_SOLVE_HH
#define D6_NEWTONIAN_SOLVE_HH


#include "utils/alg.hh"

#include <vector>

namespace d6 {

class DynParticles ;
struct Phase ;

class RigidBody ;
struct RigidBodyData ;

struct Config ;
class Stats ;

struct PhaseStepData ;

struct NewtonianSolver
{

static void solveIncompressibility(
        const Config &c, const Scalar dt, const PhaseStepData& stepData,
        std::vector<RigidBodyData> &rbData,
        DynVec &u, Phase& phase, Stats &simuStats ) ;

} ;

}

#endif
