#!/bin/bash
#
#SBATCH -J MovingCircleRbf
#SBATCH -o /home/hpc/pr63so/ga39puw2/jobs/logs/MovingCircleRbf.%j.%N.out
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
output="$(pwd)/MovingCircleRbf"
mkdir -p $output
mkdir -p $output/Precice
mkdir -p $output/Fluid
mkdir -p $output/Petsc
cp -f $conf/Precice/MovingCircle-NoIbMapping.xml $output/Precice/MovingCircle-NoIbMapping.xml
cp -f $conf/Fluid/MovingCircleRbf.xml $output/Fluid/MovingCircleRbf.xml
cp -f $conf/Structure/UnsteadyMovingCircle.xml $output/Structure/UnsteadyMovingCircle.xml
cp -f $conf/Petsc/Basic.conf $output/Petsc/Basic.conf
cd $output/Precice
mpiexec -n 1 binprecice server Fluid MovingCircle-NoIbMapping &
mpiexec -n 1 binprecice server Structure MovingCircleRbf-NoIbMapping.xml &
cd $output

mpiexec -n 1 $bin/Structure \
  -s Structure/UnsteadyMovingCircle.xml &

mpiexec -n 16 $bin/Fluid \
  -o . \
  -e Petsc/Basic.conf \
  -s Fluid/MovingCircleRbf.xml
