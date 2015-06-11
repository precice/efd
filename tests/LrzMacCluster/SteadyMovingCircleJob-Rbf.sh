#!/bin/bash
#
#SBATCH -J SteadyMovingCircle-Rbf
#SBATCH -o /home/hpc/pr63so/ga39puw2/jobs/logs/SteadyMovingCircle-Rbf.%j.%N.out
#SBATCH -D  /home/hpc/pr63so/ga39puw2/jobs/output
#SBATCH --mail-type=end
#SBATCH --mail-user=SlavaMikerov@gmail.com
#SBATCH --time=24:00:00
#SBATCH --get-user-env
#SBATCH --partition=bdz
#SBATCH --ntasks=19

DIR=$HOME/Project
bin=$DIR/.install/Gcc/Release/bin
conf=$DIR/tests/LrzMacCluster

cwd="$(pwd)"
output="$(pwd)/SteadyMovingCircle-Rbf"
mkdir -p $output
mkdir -p $output/Precice
mkdir -p $output/Fluid
mkdir -p $output/Structure
mkdir -p $output/Petsc
cp -f $conf/Precice/MovingCircle-NoIbMapping.xml $output/Precice/MovingCircle-NoIbMapping.xml
cp -f $conf/Fluid/SteadyMovingCircle-Rbf.xml $output/Fluid/SteadyMovingCircle-Rbf.xml
cp -f $conf/Structure/SteadyMovingCircle-Rbf.xml $output/Structure/SteadyMovingCircle-Rbf.xml
cp -f $conf/Petsc/Basic.conf $output/Petsc/Basic.conf
cd $output/Precice
mpiexec -n 1 binprecice server Fluid MovingCircle-NoIbMapping.xml &
mpiexec -n 1 binprecice server Structure MovingCircle-NoIbMapping.xml &
cd $output

mpiexec -n 1 $bin/Structure \
  -s Structure/SteadyMovingCircle-Rbf.xml &

mpiexec -n 16 $bin/Fluid \
  -o . \
  -e Petsc/Basic.conf \
  -s Fluid/SteadyMovingCircle-Rbf.xml
