#ifndef D6_DIPHASIC_FRICTION_HH
#define D6_DIPHASIC_FRICTION_HH

#include "DiphasicPrimal.hh"

#include "solve/FrictionSolver.hh"

#include <bogus/Core/Utils/Timer.hpp>

namespace d6 {


struct DiphasicFrictionSolver {

	struct Options
	{
		enum Algorithm {
			PG,
			GS,
			ADMM
		};

		Algorithm algorithm ;

		bool reduced   ;
		bool direct    ;
		bool useCadoux ;

		unsigned maxIterations ;	  //Inner
		unsigned maxOuterIterations ; //For Cadoux algorithm

		int projectedGradientVariant ;

		bool useInfinityNorm ;
		Scalar tolerance ;

		Options() ;
	};

	typedef Eigen::SparseMatrix< Scalar > ESM ;

	explicit DiphasicFrictionSolver( const DiphasicPrimalData& data)
		:m_data(data)
	{}


	Scalar solve( const Options& options, DynVec &x, DynVec &lambda,
				  FrictionSolver::Stats& stats	) const ;
	// x = u w p
	// M (u w) = l + (G H)' lambda
	// gamma = Gu + Hw + k
	// (lambda, gamma) \in DPmu

	Scalar solve( const Options& options, const ESM &M, const DiphasicPrimalData::MInvType& M_inv,
				  DynVec &x, DynVec &lambda, FrictionSolver::Stats& stats
				) const ;




private:

	Scalar solveRed(const Options &options, const Scalar pen,
		DynVec &x, DynVec &lambda, FrictionSolver::Stats &stats ) const ;

	Scalar solveStokes( const Options& options, const DiphasicPrimalData::MInvType& M_inv,
				  DynVec &x, DynVec &lambda, FrictionSolver::Stats& stats
				) const ;


	const DiphasicPrimalData& m_data ;

};


} //d6

#endif
