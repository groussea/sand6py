#ifndef D6_UNITS_HH
#define D6_UNITS_HH

#include "scalar.hh"

#include <cmath>

namespace d6 {

struct Units {

	enum Unit {
		None,
		Length,
		Volume,
		VolumicMass,
		Acceleration,
		Velocity,
		Time,
		Frequency,
		Viscosity,
		Stress
	} ;

	Scalar L ;  //!< Typical length (m)
	Scalar G ;  //!< Typical acceleration (m.s^{-2})
	Scalar R ;  //!< Typical density (kg.m^{-3})

	Scalar U ;  //!< Typical velocity (m.s^{-1})
	Scalar T ;  //!< Typical time (s)
	Scalar P ;  //!< Typical stress (Pa)
	Scalar M ;  //!< Typical viscosity (Pa.s)

	Units() ;

	void setTypical( Scalar length, Scalar acc, Scalar volMass ) ;

	Scalar fromSI( Unit u ) const ;
	Scalar   toSI( Unit u ) const ;

};


}

#endif
