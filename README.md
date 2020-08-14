# enriched generalized Lotka-Volterra equations
This code includes source and postprocessing code for the reduced and enriched GLV  models.
The parameters of the discrepancy model are calibrated using QUESO.

Dependencies: QUESO (https://github.com/libqueso/queso), Eigen

My Makefile is included for reference. A bin and build file should be created locally.

The steps to run the code are as follows:
```
python3 inputs/interaction-matrix.py
./bin/gen_data
rm -r outputData
./bin/glv_plus inputs/mhInput.inp
```
The first creates the interaction matrix $A$, the random vector $b$, and outputs the files info.txt and matrix.txt.

The second generates data based on the information given above.

The third removes existing 'outputData' directory (of course, alternatively you
can just rename this to save results). This directory is created during the
inverse problem.  'mhInput.inp' is an input file for QUESO.

The fourth runs the inverse problem according to the queso inputs in the mhInput.inp.

To plot the time series results:
```
cd postprocessing
./post_proc.sh
python3 compute-percents.py
python3 percent-time-series.py
```

You can also run 
```
./bin/gen_data
```
to test the forward model alone.

Notes:  
You can also run the above steps with a fixed test case by
```test_run_data_and_ip.sh
```
