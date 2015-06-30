#!/bin/bash
#
#SBATCH -J UnsteadyMovingCircle-Precice
#SBATCH -o /home/hpc/pr63so/ga39puw2/jobs/logs/UnsteadyMovingCircle-Precice.%j.%N.out
#SBATCH -D  /home/hpc/pr63so/ga39puw2/jobs/output
#SBATCH --mail-type=end
#SBATCH --mail-user=SlavaMikerov@gmail.com
#SBATCH --time=24:00:00
#SBATCH --get-user-env
#SBATCH --partition=bdz
#SBATCH --ntasks=18

DIR=$HOME/Project
bin=$DIR/.install/Gcc/Release/bin
conf=$DIR/tests/LrzMacCluster

cwd="$(pwd)"
output="$(pwd)/UnsteadyMovingCircle-Precice"
mkdir -p $output
mkdir -p $output/Precice
mkdir -p $output/Fluid
mkdir -p $output/Structure
mkdir -p $output/Petsc
cp -f $conf/Precice/MovingCircle-IbMapping.xml $output/Precice/MovingCircle-IbMapping.xml
cp -f $conf/Precice/CouplingModeDirectForcingAction.py $output/Precice/CouplingModeDirectForcingAction.py
cp -f $conf/Fluid/UnsteadyMovingCircle-Precice.xml $output/Fluid/UnsteadyMovingCircle-Precice.xml
cp -f $conf/Structure/UnsteadyMovingCircle-Precice.xml $output/Structure/UnsteadyMovingCircle-Precice.xml
cp -f $conf/Petsc/Basic.conf $output/Petsc/Basic.conf
cd $output/Precice
mpiexec -n 1 binprecice server Fluid MovingCircle-IbMapping.xml &
cd $output

mpiexec -n 1 $bin/Structure \
  -s Structure/UnsteadyMovingCircle-Precice.xml &

mpiexec -n 16 $bin/Fluid \
  -o . \
  -e Petsc/Basic.conf \
  -s Fluid/UnsteadyMovingCircle-Precice.xml
