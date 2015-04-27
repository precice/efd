#! /usr/env bash
#
module unload python
module unload pythonLib
export PYTHON_BASE=$HOME/install/python
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
export BOOST_ROOT=~/builds/install
export PRECICE_BOOST_ROOT=~/builds/install
export PRECICE_BOOST_LIB_PATH=~/builds/install/lib

module unload mpi.ibm
module unload mpi.intel
module load mpi.intel/4.1
export MPI_HOME=$MPI_BASE

module load hdf5/mpi/1.8
export HDF5_ROOT=$HDF5_BASE

module load petsc/3.5
export PETSC_DIR=/lrz/sys/libraries/petsc/3.5.2/real_mpi.intel_140_opt
export PETSC_ARCH=linux-gnu-intel

export XDMF_ROOT=~/builds/install
export EIGEN_INCLUDE_DIRS=~/builds/install/include
export PATH=~/builds/install/bin:~/install/python/bin:~/Precice/build/release:$PATH
export PRECICE_DIR=~/Precice
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/builds/install/lib:$HOME/install/python/lib
