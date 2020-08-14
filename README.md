# enriched generalized Lotka-Volterra equations
This code includes source and postprocessing code for the reduced and enriched GLV  models.
The parameters of the discrepancy model are calibrated using QUESO.

Dependencies: QUESO (https://github.com/libqueso/queso), Eigen

My Makefile is included for reference. A bin and build file should be created locally.

The steps to run the code are as follows:
```
rm -r outputData
./bin/zika_ip inputs/mhInput.inp
```
Remove or move existing 'outputData' directory. This directory is created during the inverse problem. 
'mhInput.inp' is an input file for QUESO.

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
The standard part of the code---the inverse problem---can be run with the following commands:

1) python inputs/interaction-matrix.py
2) ./bin/gen_data
3) ./bin/glv_plus inputs/mhInput.inp

The first creates the interaction matrix $A$, the random vector $b$, and outputs the files info.txt and matrix.txt.

The second generates data based on the information given above.

The third runs the inverse problem according to the queso inputs in the mhInput.inp.


Notes:  
You can also run the above steps with a fixed test case by
```test_run_data_and_ip.sh'''

<!--- To cite this code, you can use: 
```
@misc{morrison2020zikacode,  
  doi = {10.5281/ZENODO.3666845},  
  url = {https://zenodo.org/record/3666845},  
  author = {Morrison, Rebecca E.},  
  title = {rebeccaem/zika: Initial release},  
  publisher = {Zenodo},  
  version = {0.1.0},   
  year = {2020}  
 }
 ```
---!>
