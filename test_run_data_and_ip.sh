#!/bin/bash
rm -r outputData
python3 inputs/test-matrix.py
./bin/gen_data
./bin/glv_plus inputs/mhInput.inp
cd postprocessing
./post_proc.sh
