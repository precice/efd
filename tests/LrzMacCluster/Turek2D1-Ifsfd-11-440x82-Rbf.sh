#!/bin/bash
#
#SBATCH -J SeqTurek2D
#SBATCH -o /home/hpc/pr63so/ga39puw2/jobs/logs/SeqTurek2D.%j.%N.out
#SBATCH -D  /home/hpc/pr63so/ga39puw2/jobs/output
#SBATCH --mail-type=end
#SBATCH --mail-user=SlavaMikerov@gmail.com
#SBATCH --time=24:00:00
#SBATCH --get-user-env
#SBATCH --partition=bdz
#SBATCH --ntasks=1

DIR=$HOME/Project
bin=$DIR/.install/Gcc/Release/bin
conf=$DIR/tests/LrzMacCluster

source $conf/environment-configuration.sh

mpiexec -n 1 $bin/Fluid \
  -o Turek2D1-Ifsfd-11-440x82-Rbf \
  -e $conf/Petsc/Basic.conf \
  -s $conf/Fluid/Turek2D1-Ifsfd-11-440x82-Rbf.xml
