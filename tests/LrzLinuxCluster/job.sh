#!/bin/bash
#
#SBATCH -J myjob
#SBATCH -o /home/hpc/pr63so/ga39puw2/jobs/logs/myjob.%j.%N.out
#SBATCH -D  /home/hpc/pr63so/ga39puw2/jobs/output
#SBATCH --mail-type=end
#SBATCH --mail-user=SlavaMikerov@gmail.com
#SBATCH --time=05:00:00
#SBATCH --get-user-env
#SBATCH --clusters=mpp1
#SBATCH --ntasks=1

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

source $DIR/environment-configuration.sh

bin=$DIR/../../.install/Release/bin
conf=$DIR

mpiexec -n 1 $bin/Fluid \
  -o Turek2D1-Ifsfd-11-220x205-11 \
  -e $conf/Petsc/Basic.conf \
  -s $conf/Fluid/Turek2D1-Ifsfd-11-220x205-11.xml \
  -p $conf/Precice/Turek2D.xml

mpiexec -n 1 $bin/Fluid \
  -o Turek2D2-Ifsfd-11-220x205-11 \
  -e $conf/Petsc/Basic.conf \
  -s $conf/Fluid/Turek2D2-Ifsfd-11-220x205-11.xml \
  -p $conf/Precice/Turek2D.xml
