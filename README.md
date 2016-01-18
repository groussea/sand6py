
# About

This software contains a minimalist C++ implementation of the algorithms described in the document
*A Semi-Implicit Material Point Method for the Continuum Simulation of Granular Materials*.

## Contents

This archive is organized as follow:

  - `src`    Source code of the shared library
    - `utils`  Generic utilities (String and file manipulation, etc)
    - `geo`    Geometric tools (Meshes, fields, tensors)
    - `solve`  Tools for solving LCP and DCFP
	- `simu`   Proper simulation 
	- `visu`   Tools for visualizing simulation results
	- `gl`     OpenGL drawing code
  - `apps`   Source code of individual applications
  - `scenes` Sample configuration files
  - `tests`  Unit tests
  - `cmake`  CMake modules
  - `vendor` Copies of external header-only libraries, [Eigen] [1] and [So-bogus] [2] 


## Relevant files
At this stage, the documentation is still very scarce; however, sections
that are the most relevant to the article have been more densely annotated.

The reader interested in the implementation of the main simulation loop should start with the method `Simulation::step()` defined in the file `src/simu/Simu.cc`. 
Methods for reading data from particles, splitting, merging and moving them are can be found in `simu/DynParticles.cc`.
The code for computing the end-of-step velocities and stresses lies mostly inside
the `PhaseSolver::solve()` function from the file `simu/PhaseSolver.cc`.
Bindings for the external DFCP solver are in `solve/Primal.cc`.

# Usage

## Building the software

### Requirements
Building this software requires CMake and a C++11 compiler, and has
only be tested on GNU/Linux with recent versions of `g++` and `clang++`.

Other dependencies are:

  - [Eigen][1] and [So-bogus][2], which are already included in this archive
  - [boost\_serialization][3]
  - [libQGLViewer][4] (optional, required for the compilation of the OpenGL viewer)

### Compiling

The standard compilation procedure consists in typing the following
commands from the archive's root directory

	> mkdir build
	> cd build
	> cmake ..
	> make
	> ./tests/testd6

## Applications

Successful compilation should generate the following binaries in the `apps` subfolder of the `build` directory:

  - `d6` simulation tool
  - `d62vtk` Transform `d6` output to VTK files that can be read with e.g. `paraview`
  - `d6gl` OpenGL viewer and rendering utility (requires libQGLViewer)
  - `d6_solvePrimal` offline DCFP solver, for benchmarking purposes
  - `d6_analyze` post-processing tool

Usage information for these applications can be obtained with the `-?` flag.

In a typical workflow, the simulator would be run by typing e.g.

	> ./apps/d6 -i ../scenes/collapse.conf

This creates an 'out' directory, with one subdirectory per simulation frame.
If the `d6gl` application was compiled, it can be used to visualize the results:

	> ./apps/d6gl
	> ./apps/d6gl -s 16 -g .2    #Render grain samples instead of particules

Pressing 'H' inside the viewer triggers an help window listing libQGLViewer mouse and keyboard commands.
A useful yet unintuitive shortcut is the 'double-middle-click' which centers the scene.

## Configuration fields

Here is a list of the various configuration options that can be passed to the `d6` application, with their default values. Configuration files corresponding to the simulations presented in the article are also provided in the `scenes` directory, and can be used as e.g. `./apps/d6 -i ../scenes/scene_file.conf`.

### Simulation size and resolution

  - `fps`=*240 s^{-1}*  Number of frames to be generated per SI seconds
  - `substeps`=*1*      Number of simulation sub-steps per frame. `0` means adaptive.
  - `nFrames`=*1*       Number of frames to generate
  - `box`=*(1,1,1) m*   Grid domain
  - `res`=*(10,10,10)*  Number of grid cells for each dimension
  - `nSamples`=*2*      Number of particles to generate per grid cell per dimension
  - `randomize`=*0*     Amount of random perturbation to add to initial particle positions
  - `boundary`=*cuve*   Boundary conditions ( stick, slip, ... ) preset 
  - `scenario`=*bed*    Scene configuration preset

### Physical parameters

  - `volMass`=*1.5e3 kg.m^{-3}*     Volumetric mass density of the grains
  - `viscosity`=*1.e-3 Pa s*        Newtonian dynamic viscosity 
  - `gravity`=*(0,0,-9.81) m^{-2}*  Acceleration due to gravity
  - `phiMax` = *1*				    Maximum volume fraction of grains
  - `mu`=*0*                     Grains friction coefficient
  - `delta_mu`=*0*               Dynamic velocity increase for mu(I) rheology
  - `I0`=*0.4*                   Inertial number for mu(I) rheology
  - `grainDiameter`=*1.e-3 m*    Grain diameter for mu(I) rheology
  - `muRigid`=*0*				 Grain/body friction coefficient
  - `cohesion`=*0 Pa*		 Grain cohesion
  - `cohesion_decay`=*0*     Cohesion decay rate in shear flows
  - `anisotropy`=*0*         Amount of friction anisotropy
  - `elongation`=*0*         Tendency of particles to align with flow
  - `brownian`=*0*	         Tendency of particles to return to isotropic orientation
  - `initialOri`=*(1/3,1/3,1/3)* Initial eigenvalues of orientation tensor in grid-aligned basis

### Miscellaneous

  - `enforceMaxFrac`=*false*  Perform geometric volume correction projection
  - `output`=*true*       Output simulation states (field and particle data)
  - `dumpPrimalData`=*0*  If non-zero, dump every *n*th primal problems
  


# License

This software is distributed under the terms of the [GNU General Public License Version 3][5].

  [1]: http://eigen.tuxfamily.org     "Eigen, template library for linear algebra"
  [2]: http://gdaviet.fr/code/bogus   "So-bogus, Coulomb friction solver"
  [3]: http://www.boost.org/doc/libs/release/libs/serialization/ "Boost serialization library"
  [4]: http://libqglviewer.com        "Qt-base OpenGL viewer framework"
  [5]: http://www.gnu.org/licenses/gpl-3.0.en.html "GNU General Public License Version 3"
