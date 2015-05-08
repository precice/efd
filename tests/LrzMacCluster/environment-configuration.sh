#! /usr/env bash

module unload python
module unload pythonLib
export PYTHON_BASE=$HOME/install/icc-15.0/python-2.7.9
export PRECICE_PYTHON_INC_PATH=$PYTHON_BASE/include/python2.7
export PRECICE_NUMPY_INC_PATH=$PYTHON_BASE/lib/python2.7/site-packages/numpy/core/include/numpy
export PRECICE_PYTHON_LIB_PATH=$PYTHON_BASE/lib
export PRECICE_SOCKET_LIB=pthread

module load git
module unload gcc
module load gcc/4.9
module unload ccomp/intel
module load ccomp/intel/15.0
module unload cmake

module unload boost
export BOOST_ROOT=~/install/icc-15.0/boost-1.57
export PRECICE_BOOST_ROOT=$BOOST_ROOT
export PRECICE_BOOST_LIB_PATH=$BOOST_ROOT/lib

module unload mpi.ibm
module unload mpi.intel
module load mpi.intel/4.1
export MPI_HOME=$MPI_BASE

module load hdf5/mpi/1.8
export HDF5_ROOT=$HDF5_BASE

module load petsc/3.5
export PETSC_DIR=/lrz/sys/libraries/petsc/3.5.2/real_mpi.intel_140_opt
export PETSC_ARCH=linux-gnu-intel

export EIGEN_INCLUDE_DIRS=~/install/icc-15.0/eigen-3.2.4/include
export PATH=$HOME/install/icc-15.0/Uni/bin:$MPI_BASE/intel64/bin:$PYTHON_BASE/bin:~/install/icc-15.0/cmake-3.2.2/bin:~/Precice/build/release:$PATH
export PRECICE_DIR=~/Precice
export LD_LIBRARY_PATH=$MPI_BASE/intel64/lib:$PYTHON_BASE/lib:$BOOST_ROOT/lib:$LD_LIBRARY_PATH
export CXX="icc"
export CC="icc"
export LDFLAGS="-static-intel -static-libgcc"
