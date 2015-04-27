#!/bin/bash
#
#SBATCH -J myjob
#SBATCH -o /home/hpc/pr63so/ga39puw2/jobs/logs/myjob.%j.%N.out
#SBATCH -D  /home/hpc/pr63so/ga39puw2/jobs/output
#SBATCH --mail-type=end
#SBATCH --mail-user=SlavaMikerov@gmail.com
#SBATCH --time=24:00:00
#SBATCH --get-user-env
#SBATCH --partition=snb
#SBATCH --ntasks=1

DIR=$HOME/Project
bin=$DIR/.install/Release/bin
conf=$DIR/tests/LrzMacCluster

source $conf/environment-configuration.sh

mpiexec -n 1 $bin/Fluid \
  -o Turek2D1-Ifsfd-11-440x82-11 \
  -e $conf/Petsc/Basic.conf \
  -s $conf/Fluid/Turek2D1-Ifsfd-11-440x82-11.xml \
  -p $conf/Precice/Turek2DNoIb-440x82.xml

mpiexec -n 1 $bin/Fluid \
  -o Turek2D2-Ifsfd-11-440x82-11 \
  -e $conf/Petsc/Basic.conf \
  -s $conf/Fluid/Turek2D2-Ifsfd-11-440x82-11.xml \
  -p $conf/Precice/Turek2DNoIb-440x82.xml
