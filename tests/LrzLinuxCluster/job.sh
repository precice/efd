#!/bin/bash
#
#SBATCH -J myjob
#SBATCH -o /home/hpc/pr63so/ga39puw2/jobs/logs/myjob.%j.%N.out
#SBATCH -D  /home/hpc/pr63so/ga39puw2/jobs/output
#SBATCH --mail-type=end
#SBATCH --mail-user=SlavaMikerov@gmail.com
#SBATCH --time=01:00:00
#SBATCH --get-user-env
#SBATCH --clusters=mpp1
#SBATCH --ntasks=25

source ./modules.sh

bin=~/Project/.install/Release/bin
conf=~/Project/tests/LrzLinuxCluster

mpiexec -n 25 $bin/Fluid \
  -o Turek2D1-Ifsfd-11-220x205 \
  -e $conf/Petsc/Basic.conf \
  -s $conf/Fluid/Turek2D1-Ifsfd-11-220x205.xml \
  -p $conf/Precice/Turek2D.xml
