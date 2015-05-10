#!/usr/bin/env python
# encoding: utf-8

import os
import re
import subprocess
import shutil
import sys

# Add and parse arguments {{{
# -----------------------------------------------------------------------------
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument('-f', nargs='?', default=".*")
parser.add_argument('-p', nargs='?', default="1")
parser.add_argument('-d', nargs=1)

args, unknown_args = parser.parse_known_args()
# -----------------------------------------------------------------------------
# }}}

template = (""
"#!/bin/bash\n"
"#\n"
"#SBATCH -J {0}\n"
"#SBATCH -o /home/hpc/pr63so/ga39puw2/jobs/logs/{0}.%j.%N.out\n"
"#SBATCH -D  /home/hpc/pr63so/ga39puw2/jobs/output\n"
"#SBATCH --mail-type=end\n"
"#SBATCH --mail-user=SlavaMikerov@gmail.com\n"
"#SBATCH --time=24:00:00\n"
"#SBATCH --get-user-env\n"
"#SBATCH --partition=bdz\n"
"#SBATCH --ntasks={1}\n"
"\n"
"DIR=$HOME/Project\n"
"bin=$DIR/.install/Gcc/Release/bin\n"
"conf=$DIR/tests/LrzMacCluster\n"
"\n"
"source $conf/environment-configuration.sh\n"
"\n"
"mpiexec -n {1} $bin/Fluid \\\n"
"  -o {0} \\\n"
"  -e $conf/Petsc/Basic.conf \\\n"
"  -s $conf/Fluid/{0}.xml")

for root, dirs, files in os.walk(args.d[0]):
  for file_name in files:
    if not re.match(args.f, file_name):
      continue
    basic_file_name = os.path.splitext(file_name)[0]
    sh_file_name = basic_file_name + ".sh"
    if os.path.exists(sh_file_name):
      os.remove(sh_file_name)
    fo = open(sh_file_name, "wb+")
    fo.write(template.format(basic_file_name, args.p));
    fo.close()
    subprocess.call("sbatch {0} {1}".format(" ".join(unknown_args), sh_file_name), shell=True)
