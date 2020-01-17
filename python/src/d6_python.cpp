#include <pybind11/pybind11.h>
#include "utils/Config.hh"
#include "mono/MonoSimu.hh"

#include "utils/Log.hh"
#include "utils/File.hh"

#include <cstring>

#include "utils/string.hh"

#include "visu/Offline.hh"
#include "visu/VTKParticlesWriter.hh"
#include "visu/VTKFieldWriter.hh"

#include "mono/Phase.hh"

#include "geo/FieldBase.impl.hh"
#include "geo/ScalarField.hh"
#include "geo/TensorField.hh"

#include <iostream>




void dump_frame( const d6::Offline& offline, bool particles,
                 const char* base_dir, unsigned frame )
{

	if(particles) {
		d6::VTKParticlesWriter particlesWriter( base_dir, offline.particles() ) ;
		particlesWriter.startFile( "particles", frame ) ;
		particlesWriter.dump_all() ;
	}


		d6::VTKFieldWriter<d6::PrimalShape> fieldWriter( base_dir, offline.meshes().primal() ) ;
	//	fieldWriter.setMode( d6::VTKWriter::Ascii );
		fieldWriter.startFile( "primal-fields", frame ) ;
		fieldWriter.dump(    "phi", offline.grains().fraction ) ;
		fieldWriter.dump(      "u", offline.grains().velocity ) ;
		fieldWriter.dump(  "d_phi", offline.grains().grad_phi ) ;
		if( offline.config().exportAllFields )
		{
			fieldWriter.dump( "forces", offline.grains().fcontact ) ;

			d6::PrimalTensorField tau = offline.grains().stresses.interpolate<d6::PrimalShape>(
			            offline.meshes().primal()) ;
			fieldWriter.dump("stresses", tau ) ;
		}
	

	if( offline.config().exportAllFields )
	{
 #ifdef D6_UNSTRUCTURED_DUAL
 		d6::VTKParticlesWriter fieldWriter( base_dir, offline.particles() ) ;
 #else
		 d6::VTKFieldWriter<d6::DualShape> fieldWriter( base_dir, offline.meshes().primal() ) ;
 #endif

		 fieldWriter.startFile( "dual-fields", frame ) ;

		d6::DualScalarField p   = offline.grains().stresses.trace() ;
		d6::DualScalarField dh  = offline.grains().sym_grad.trace() ;
		d6::DualScalarField taun= d6::DualTensorField( offline.grains().stresses.deviatoricPart() ).norm() ;
		fieldWriter.dump(    "p", p) ;

		fieldWriter.dump(   "dh", dh ) ;

		fieldWriter.dump( "taun", taun ) ;

			fieldWriter.dump( "lambda", offline.grains().stresses ) ;
	}
}


int d62vtk(const char* base_dir = "out" ,bool allF = false , uint frame = 0, bool particles =false) {


    d6::Offline offline(base_dir) ;

	unsigned cur_frame = frame ;

	do {
		if(! offline.load_frame( cur_frame ) )
			return allF?0:1 ;
		dump_frame( offline, particles, base_dir, cur_frame++ ) ;
	} while( allF ) ;

	return (frame == cur_frame) ? 1 : 0;
    
    };


namespace py = pybind11;

PYBIND11_MODULE(d6_python, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: d6_python

        .. autosummary::

    )pbdoc";


    m.def("d6run",[](const char* base_dir,const char* config_file) {	

    d6::Config config ;
    config.from_file(config_file);


	//write base_dir on config
	config.from_string("base_dir", base_dir);
	// Save copy of final configuration and convert to interal units
	d6::FileInfo outDir ( base_dir ) ;
	if( !outDir.exists() ) outDir.makeDir() ;
	config.dump( outDir.filePath("config"), "config_out" );
	config.internalize();

	d6::Log::Debug() << "Typical length = " << config.units().toSI(d6::Units::Length) << " m"<< std::endl ;
	d6::Log::Debug() << "Typical velocity = " << config.units().toSI(d6::Units::Velocity) << " m.s^-1" << std::endl ;
	d6::Log::Debug() << "Typical pressure = " << config.units().toSI(d6::Units::Stress) << " Pa" << std::endl ;
	d6::Log::Debug() << "1/Re = " << config.viscosity << std::endl ;

	// Run simulation
	d6::MonoSimu( config, base_dir ).run() ;  
    return 0 ;
        } );

    m.def("d62vtk", &d62vtk , py::arg("base_dir") = py::cast("out") ,  
	py::arg("allF") = false, py::arg("frame") = 0, py::arg("particles") = false  
	);


#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
