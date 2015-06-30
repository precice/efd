#!/bin/bash
#
#SBATCH -J Unsteady100MovingCircle-Rbf
#SBATCH -o /home/hpc/pr63so/ga39puw2/jobs/logs/Unsteady100MovingCircle-Rbf.%j.%N.out
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
output="$(pwd)/Unsteady100MovingCircle-Rbf"
mkdir -p $output
mkdir -p $output/Precice
mkdir -p $output/Fluid
mkdir -p $output/Structure
mkdir -p $output/Petsc
cp -f $conf/Precice/MovingCircle-NoIbMapping.xml $output/Precice/MovingCircle-NoIbMapping.xml
cp -f $conf/Fluid/Unsteady100MovingCircle-Rbf.xml $output/Fluid/Unsteady100MovingCircle-Rbf.xml
cp -f $conf/Structure/Unsteady100MovingCircle-Rbf.xml $output/Structure/Unsteady100MovingCircle-Rbf.xml
cp -f $conf/Petsc/Basic.conf $output/Petsc/Basic.conf
cd $output/Precice
mpiexec -n 1 binprecice server Fluid MovingCircle-NoIbMapping.xml &
mpiexec -n 1 binprecice server Structure MovingCircle-NoIbMapping.xml &
cd $output

mpiexec -n 1 $bin/Structure \
  -s Structure/Unsteady100MovingCircle-Rbf.xml &

mpiexec -n 16 $bin/Fluid \
  -o . \
  -e Petsc/Basic.conf \
  -s Fluid/Unsteady100MovingCircle-Rbf.xml
