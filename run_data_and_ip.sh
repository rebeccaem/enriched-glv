#!/bin/bash
rm -r outputData
python3 inputs/interaction-matrix.py
./bin/gen_data
./bin/glv_plus inputs/mhInput.inp
cd postprocessing
./post_proc.sh
#python compute-gammas.py
#python append-point-gammas.py
