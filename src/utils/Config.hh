#ifndef D6_CONFIG_HH
#define D6_CONFIG_HH

#include "utils/units.hh"
#include "utils/alg.hh"

#include <string>

namespace d6 {


#define EXPAND_CONFIG \
	CONFIG_FIELD( fps				, Scalar		,	Units::Frequency	) \
	CONFIG_FIELD( substeps			, unsigned		,	Units::None			) \
	CONFIG_FIELD( nFrames			, unsigned		,	Units::None			) \
	\
	CONFIG_FIELD( box				, Vec			,	Units::Length		) \
	CONFIG_FIELD( res				, Vec3i			,	Units::None			) \
	CONFIG_FIELD( nSamples			, unsigned		,	Units::None			) \
	\
	CONFIG_FIELD( volMass			, Scalar		,	Units::VolumicMass	) \
	CONFIG_FIELD( viscosity			, Scalar		,	Units::Viscosity	) \
	CONFIG_FIELD( gravity			, Vec			,	Units::Acceleration	) \
	CONFIG_FIELD( phiMax			, Scalar		,	Units::None			) \
	CONFIG_FIELD( mu				, Scalar		,	Units::None			) \
	CONFIG_FIELD( cohesion			, Scalar		,	Units::Stress		) \
	CONFIG_FIELD( cohesion_decay	, Scalar		,	Units::None			) \
	\
	CONFIG_FIELD( scenario			, std::string	,	Units::None			) \
	CONFIG_FIELD( boundary			, std::string	,	Units::None			) \
	CONFIG_FIELD( output			, bool			,	Units::None			) \


struct Config
{
	Config( ) ;

	bool from_string( const std::string &key, const std::string &value ) ;
	bool from_file(const std::string& file_name) ;

	bool dump( const std::string& file_name ) const ;

	//! Transform all the parameters from SI to internal units
	void internalize() ;

	const Units& units() const {
		return m_units ;
	}

	Scalar inv_dt() const {
		return ( fps * substeps ) ;
	}
	Scalar dt() const {
		return 1./ inv_dt() ;
	}

#define CONFIG_FIELD( name, type, u ) \
	type name ;
	EXPAND_CONFIG
#undef CONFIG_FIELD

private:
	bool from_string( const std::string &key, std::istringstream &value ) ;

	Units m_units ;

} ;

} //ns hyb2d

#endif

