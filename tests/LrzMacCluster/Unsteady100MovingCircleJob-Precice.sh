#!/bin/bash
#
#SBATCH -J Unsteady100MovingCircle-Precice
#SBATCH -o /home/hpc/pr63so/ga39puw2/jobs/logs/Unsteady100MovingCircle-Precice.%j.%N.out
#SBATCH -D  /home/hpc/pr63so/ga39puw2/jobs/output
#SBATCH --mail-type=end
#SBATCH --mail-user=SlavaMikerov@gmail.com
#SBATCH --time=24:00:00
#SBATCH --get-user-env
#SBATCH --partition=bdz
#SBATCH --ntasks=3

DIR=$HOME/Project
bin=$DIR/.install/Gcc/Release/bin
conf=$DIR/tests/LrzMacCluster

cwd="$(pwd)"
output="$(pwd)/Unsteady100MovingCircle-Precice"
mkdir -p $output
mkdir -p $output/Precice
mkdir -p $output/Fluid
mkdir -p $output/Structure
mkdir -p $output/Petsc
cp -f $conf/Precice/MovingCircle-IbMapping.xml $output/Precice/MovingCircle-IbMapping.xml
cp -f $conf/Precice/CouplingModeDirectForcingAction.py $output/Precice/CouplingModeDirectForcingAction.py
cp -f $conf/Fluid/Unsteady100MovingCircle-Precice.xml $output/Fluid/Unsteady100MovingCircle-Precice.xml
cp -f $conf/Structure/Unsteady100MovingCircle-Precice.xml $output/Structure/Unsteady100MovingCircle-Precice.xml
cp -f $conf/Petsc/Basic.conf $output/Petsc/Basic.conf
cd $output/Precice
# mpiexec -n 1 binprecice server Fluid MovingCircle-IbMapping.xml &
cd $output

mpiexec -n 1 $bin/Structure \
  -s Structure/Unsteady100MovingCircle-Precice.xml &

mpiexec -n 1 $bin/Fluid \
  -o . \
  -e Petsc/Basic.conf \
  -s Fluid/Unsteady100MovingCircle-Precice.xml
