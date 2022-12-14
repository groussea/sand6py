/*
 * This file is part of Sand6, a C++ continuum-based granular simulator.
 *
 * Copyright 2016 Gilles Daviet <gilles.daviet@inria.fr> (Inria - Universit√© Grenoble Alpes)
 *
 * Sand6 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * Sand6 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with Sand6.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef D6_PHASE_STEP_DATA_HH
#define D6_PHASE_STEP_DATA_HH

#include "PhaseFields.hh"
#include "RigidBodyData.hh"


#include "simu/ActiveIndices.hh"
#include "simu/FormBuilder.hh"
#include "simu/StressMassMatrix.hh"


namespace d6 {

struct Config ;
class RigidBody ;
class DynParticles ;
struct Phase ;

struct PhaseStepData {

	//! Active nodes: nodes under the influence of at least one particle
	Active primalNodes ;
	Active dualNodes ;

	struct Forms {
		// Linear form vectors

		DynVec phiu ;           //!< integral of fraction times velocity
		DynVec externalForces ; //!< integral of external forces
		DynArr volumes  ;	    //!< Volume aossciated to each node = int(1)

		DynArr fraction  ;		//!< integrated volume fraction occupied by the phase

		// Bilinear Form marices

		//! Lumped mass matrix, its inverse and a factorization
		typename FormMat<WD,WD>::SymType M_lumped ;
		typename FormMat<WD,WD>::SymType M_lumped_inv ;
		typename FormMat<WD,WD>::SymType M_lumped_inv_sqrt ;

		typename FormMat<WD,WD>::Type A ; //!< Mass + Visco ; Could be Symmetric when FormBuilder has sym index

		typename FormMat<SD,WD>::Type B ; //!< \phi Tau:D(u)
		typename FormMat<RD,WD>::Type J ; //!< \phi Tau:W(u)

		typename FormMat<1 ,WD>::Type C ; //!< \phi v.grad(p)


		typedef AbstractStressMassMatrix< DualShape > StressMassMatrix ;
		StressMassMatrix S ;

	} forms ;

	//! Projectors enforcing boundary conditions
	struct Projectors {
		typename FormMat<WD,WD>::SymType vel ;
		typename FormMat<SD,SD>::SymType stress ;
		typename FormMat<1,SD> ::Type    pressure;
	} proj ;


	DynArr cohesion ;   //!< interpolated cohesivity
	DynArr inertia  ;   //!< interpolated inertial number
	DynArr int_pressure ; 	//!< interpolated pressure
	DynArr DuT ; 	//!< interpolated DuT


	typename FormMat<SD,SD>::SymType Aniso ; //!< Anisotropy linear operator (N^{-1})

	PhaseStepData()
	    : m_totRbNodes( 0 )
	{}

	Index nPrimalNodes() const
	{
		return primalNodes.count() ;
	}

	Index nDualNodes() const
	{
		return dualNodes.count() ;
	}

	Index nSuppNodes() const
	{
		return m_totRbNodes ;
	}

	void compute(
	        const DynParticles& particles, const Config &config, const Scalar dt,
	        Phase &phase,
	        std::vector< RigidBody   >& rigidBodies,
	        std::vector< RBStresses > &rbStresses,
	        std::vector< RigidBodyData > &rbData  ) ;


	// Re-used for diphasic

	static void computeActiveNodes(const std::vector< bool >& activeCells,
	                               const PrimalShape &pShape , const DualShape &dShape,
	                               Active& primalNodes, Active& dualNodes
	                               ) ;

	static void computePhiAndGradPhi(const PrimalScalarField& intPhi, PrimalScalarField&phi, PrimalVectorField &grad_phi ) ;

	static void computeProjectors(const Config &config, const PrimalShape &pShape, const DualShape &dShape,
	                       const PrimalScalarField& lumped_mass,
	                       const Active& primalNodes, const Active& dualNodes, const Index nSuppNodes,
	                       Projectors& mats ) ;

	static void computeRbProjectors(const Config& config, const PrimalShape &pShape,
	                       const std::vector< RigidBodyData > &rbData, Projectors& mats ) ;

private:

	void computeActiveBodies( std::vector<RigidBody> &rigidBodies,
	                          std::vector<RBStresses> &rbStresses,
	                          std::vector< RigidBodyData > &rbData ) ;

	void computeAnisotropy(const DynVec& orientation,  const Config &config,
	                       typename FormMat<SD,SD>::SymType &Aniso ) const ;


	void assembleMatrices(const Particles& particles, const Config& c, const Scalar dt,
	                       const DualShape &dShape, const PrimalScalarField &phiInt,
	                       std::vector< RigidBodyData > &rbData ) ;


	PhaseStepData( const PhaseStepData& ) = delete ;
	PhaseStepData& operator=( const PhaseStepData& ) = delete ;

	Index m_totRbNodes ;
};


} //d6


#endif
