#!/bin/sh
bogus=https://bitbucket.org/gdaviet/so-bogus/get/master.tar.bz2
eigen=https://bitbucket.org/eigen/eigen/get/3.2.9.tar.bz2

# Build directory
mkdir pkg
cd pkg
git clone git+ssh://gdaviet@scm.gforge.inria.fr/gitroot/sand6/sand6.git

# Prepare Vendor directory
mkdir sand6/vendor
mkdir sand6/vendor/eigen3
mkdir vendor-src
cd vendor-src
mkdir bogus ; cd bogus ; wget ${bogus} ; cd ..
mkdir eigen ; cd eigen ; wget ${eigen} ; cd ..
tar -jxf eigen/*.tar.bz2 
tar -jxf bogus/*.tar.bz2
cp -r gdaviet-*/src ../sand6/vendor/bogus
cp -r gdaviet-*/LICENSE.md ../sand6/vendor/bogus/
cp -r gdaviet-*/RELEASE.md ../sand6/vendor/bogus/
rm ../sand6/vendor/bogus/CMake*
rm ../sand6/vendor/bogus/Interfaces/*.cpp
cp -r eigen-*/Eigen ../sand6/vendor/eigen3/
cp -r eigen-*/unsupported ../sand6/vendor/eigen3/
cp -r eigen-*/signature* ../sand6/vendor/eigen3/
cp -r eigen-*/COPYING.* ../sand6/vendor/eigen3/
rm -Rf gdaviet-*
rm -Rf eigen-*
# Prepare sources and make version file
cd ../sand6
git_hash=`git log -1 --format=%h`
git_branch=`git rev-parse --abbrev-ref HEAD`
build_date=`date`
rm gen.sh
rm -Rf .git
rm scenes/*.obj
echo "${git_branch}\n${git_hash}\nPackaged ${build_date}" > VERSION
# Make archive
cd ..
tar -czf sand6.tar.gz sand6
# Tests
mkdir test
cd test 
tar -xzf ../sand6.tar.gz
cd sand6
mkdir build
cd build
cmake ..
make -j4
./tests/testd6
./apps/d6

