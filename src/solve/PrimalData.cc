#include "PrimalData.hh"

#include "utils/Log.hh"

#include <bogus/Core/Block.io.hpp>
#include <bogus/Core/Block.impl.hpp>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>

#include <fstream>

namespace d6 {


template <class Archive>
void PrimalData::serialize(Archive &ar, const unsigned int version )
{
	ar & H ;
	ar & w ;

	ar & mu ;

	ar & jacobians ;
	ar & inv_inertia_matrices ;

	if( version > 0 ) {
		ar & mass_matrix_mode ;
		ar & M ;
		ar & f ;
	}
}

bool PrimalData::load(const char *file)
{
	std::ifstream ifs( file );

	if(!ifs) {
		Log::Error() << "Cannot read " << file << std::endl ;
		return false ;
	}

	boost::archive::binary_iarchive ia(ifs);

	mass_matrix_mode = Lumped ;
	ia >> (*this) ;

	return true ;
}

bool PrimalData::dump(const char *file) const
{
	std::ofstream ofs( file );

	if(!ofs) {
		Log::Error() << "Cannot write " << file << std::endl ;
		return false ;
	}

	boost::archive::binary_oarchive oa(ofs);
	oa << (*this) ;
	return true ;
}


} //d6

BOOST_CLASS_VERSION(d6::PrimalData, 1)
