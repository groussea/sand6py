#ifndef D6_STATS_HH
#define D6_STATS_HH

#include "utils/alg.hh"
#include "utils/File.hh"

namespace d6 {

#define EXPAND_STAT \
	STAT_FIELD( stepId					, unsigned		, "step"	) \
	STAT_FIELD( frameId					, unsigned		, "frame"	) \
	STAT_FIELD( delta_t					, Scalar		, "dtSI"	) \
	STAT_FIELD( nParticles				, unsigned		, "nPart"	) \
	STAT_FIELD( nNodes					, unsigned		, "totNds"  ) \
	STAT_FIELD( nActiveNodes			, unsigned		, "actNds"  ) \
	STAT_FIELD( nCouplingNodes			, unsigned		, "cplNds" ) \
	STAT_FIELD( maxVelocity				, Scalar		, "maxVel"	) \
	STAT_FIELD( frictionError			, Scalar		, "slvErr"	) \
	STAT_FIELD( frictionIterations		, unsigned		, "slvIter"	) \
	STAT_FIELD( assemblyTime			, Scalar		, "asmTime"	) \
	STAT_FIELD( linSolveTime			, Scalar		, "linTime"	) \
	STAT_FIELD( lcpSolveTime			, Scalar		, "lcpTime"	) \
	STAT_FIELD( frictionTime			, Scalar		, "slvTime"	) \
	STAT_FIELD( advectionTime			, Scalar		, "advTime"	) \
	STAT_FIELD( totalTime				, Scalar		, "totTime"	) \

class Stats {


public:
	Stats( const char* base_dir )	 ;

	void dump() ;

private:
	File 	 m_out ;

public:

#define STAT_FIELD( name, type, abv ) \
	type name ;
	EXPAND_STAT
#undef STAT_FIELD

};


} // ns d6


#endif
