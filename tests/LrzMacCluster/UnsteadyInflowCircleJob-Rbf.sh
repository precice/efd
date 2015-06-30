#!/bin/bash
#
#SBATCH -J UnsteadyInflowCircle-Rbf
#SBATCH -o /home/hpc/pr63so/ga39puw2/jobs/logs/UnsteadyInflowCircle-Rbf.%j.%N.out
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
output="$(pwd)/UnsteadyInflowCircle-Rbf"
mkdir -p $output
mkdir -p $output/Precice
mkdir -p $output/Fluid
mkdir -p $output/Structure
mkdir -p $output/Petsc
cp -f $conf/Precice/InflowCircle-NoIbMapping.xml $output/Precice/InflowCircle-NoIbMapping.xml
cp -f $conf/Fluid/UnsteadyInflowCircle-Rbf.xml $output/Fluid/UnsteadyInflowCircle-Rbf.xml
cp -f $conf/Structure/InflowCircle-Rbf.xml $output/Structure/InflowCircle-Rbf.xml
cp -f $conf/Petsc/Basic.conf $output/Petsc/Basic.conf
cd $output/Precice
# mpiexec -n 1 binprecice server Fluid InflowCircle-NoIbMapping.xml &
cd $output

mpiexec -n 1 $bin/Structure \
  -s Structure/InflowCircle-Rbf.xml &

mpiexec -n 1 $bin/Fluid \
  -o . \
  -e Petsc/Basic.conf \
  -s Fluid/UnsteadyInflowCircle-Rbf.xml
