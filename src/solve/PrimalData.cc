#include "PrimalData.hh"

#include "utils/Log.hh"

#include <bogus/Core/Block.io.hpp>
#include <bogus/Core/Block.impl.hpp>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include <fstream>

namespace d6 {


template <class Archive>
void PrimalData::serialize(Archive &ar, const unsigned int )
{
	ar & H ;
	ar & w ;

	ar & mu ;

	ar & jacobians ;
	ar & inv_inertia_matrices ;
}

bool PrimalData::load(const char *file)
{
	std::ifstream ifs( file );

	if(!ifs) {
		Log::Error() << "Cannot read " << file << std::endl ;
		return false ;
	}

	boost::archive::binary_iarchive ia(ifs);
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
