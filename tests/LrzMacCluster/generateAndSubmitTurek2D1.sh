#!/usr/bin/env bash

./generateAndSubmitJobs.py -d Fluid/ -f "Turek2D1-\w+-1x1-220x41-\w+-Precice\w*" -t 2 -m 1 --with-server -p Turek2D-IbMapping-220x41 "$@"
./generateAndSubmitJobs.py -d Fluid/ -f "Turek2D1-\w+-1x1-440x82-\w+-Precice\w*" -t 2 -m 1 --with-server -p Turek2D-IbMapping-440x82 "$@"
./generateAndSubmitJobs.py -d Fluid/ -f "Turek2D1-\w+-1x1-220x41-Rbf\w*" -t 1 -p Turek2D-NoIbMapping-220x41 "$@"
./generateAndSubmitJobs.py -d Fluid/ -f "Turek2D1-\w+-1x1-440x82-Rbf\w*" -t 1 -p Turek2D-NoIbMapping-440x82 "$@"

./generateAndSubmitJobs.py -d Fluid/ -f "Turek2D1-\w+-8x2-440x82-\w+-Precice\w*" -t 17 -m 16 --with-server -p Turek2D-IbMapping-440x82 "$@"
./generateAndSubmitJobs.py -d Fluid/ -f "Turek2D1-\w+-8x2-440x82-Rbf\w*" -t 16 -p Turek2D-NoIbMapping-440x82 "$@"

./generateAndSubmitJobs.py -d Fluid/ -f "Turek2D1-\w+-8x4-880x164-\w+-Precice\w*" -t 33 -m 32 --with-server -p Turek2D-IbMapping-880x164 "$@"
./generateAndSubmitJobs.py -d Fluid/ -f "Turek2D1-\w+-16x4-880x164-Rbf\w*" -t 64  -p Turek2D-NoIbMapping-880x164 "$@"
