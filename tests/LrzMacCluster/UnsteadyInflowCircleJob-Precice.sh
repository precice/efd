#!/bin/bash
#
#SBATCH -J UnsteadyInflowCircle-Precice
#SBATCH -o /home/hpc/pr63so/ga39puw2/jobs/logs/UnsteadyInflowCircle-Precice.%j.%N.out
#SBATCH -D  /home/hpc/pr63so/ga39puw2/jobs/output
#SBATCH --mail-type=end
#SBATCH --mail-user=SlavaMikerov@gmail.com
#SBATCH --time=24:00:00
#SBATCH --get-user-env
#SBATCH --partition=bdz
#SBATCH --ntasks=2

DIR=$HOME/Project
bin=$DIR/.install/Gcc/Release/bin
conf=$DIR/tests/LrzMacCluster

cwd="$(pwd)"
output="$(pwd)/UnsteadyInflowCircle-Precice"
mkdir -p $output
mkdir -p $output/Precice
mkdir -p $output/Fluid
mkdir -p $output/Structure
mkdir -p $output/Petsc
cp -f $conf/Precice/InflowCircle-IbMapping.xml $output/Precice/InflowCircle-IbMapping.xml
cp -f $conf/Precice/CouplingModeDirectForcingAction.py $output/Precice/CouplingModeDirectForcingAction.py
cp -f $conf/Fluid/UnsteadyInflowCircle-Precice.xml $output/Fluid/UnsteadyInflowCircle-Precice.xml
cp -f $conf/Structure/InflowCircle-Precice.xml $output/Structure/InflowCircle-Precice.xml
cp -f $conf/Petsc/Basic.conf $output/Petsc/Basic.conf
cd $output/Precice
# mpiexec -n 1 binprecice server Fluid InflowCircle-IbMapping.xml &
cd $output

mpiexec -n 1 $bin/Structure \
  -s Structure/InflowCircle-Precice.xml &

mpiexec -n 1 $bin/Fluid \
  -o . \
  -e Petsc/Basic.conf \
  -s Fluid/UnsteadyInflowCircle-Precice.xml
