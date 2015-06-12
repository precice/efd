#!/bin/bash
#
#SBATCH -J Turek3D-Rbf
#SBATCH -o /home/hpc/pr63so/ga39puw2/jobs/logs/Turek3D-Rbf.%j.%N.out
#SBATCH -D  /home/hpc/pr63so/ga39puw2/jobs/output
#SBATCH --mail-type=end
#SBATCH --mail-user=SlavaMikerov@gmail.com
#SBATCH --time=24:00:00
#SBATCH --get-user-env
#SBATCH --partition=bdz
#SBATCH --ntasks=17

DIR=$HOME/Project
bin=$DIR/.install/Gcc/Release/bin
conf=$DIR/tests/LrzMacCluster

cwd="$(pwd)"
output="$(pwd)/Turek3D-Rbf"
mkdir -p $output
mkdir -p $output/Precice
mkdir -p $output/Fluid
mkdir -p $output/Petsc
cp -f $conf/Precice/MovingCircle-NoIbMapping.xml $output/Precice/MovingCircle-NoIbMapping.xml
cp -f $conf/Fluid/Turek3D-Rbf.xml $output/Fluid/Turek3D-Rbf.xml
cp -f $conf/Petsc/Basic.conf $output/Petsc/Basic.conf
cd $output/Precice
mpiexec -n 1 binprecice server Fluid Turek3D-NoIbMapping.xml &
cd $output

mpiexec -n 16 $bin/Fluid \
  -o . \
  -e Petsc/Basic.conf \
  -s Fluid/Turek3D-Rbf.xml

