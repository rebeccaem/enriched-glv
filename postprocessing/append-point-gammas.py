import matplotlib
import numpy as np
from numpy import loadtxt
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

dataFile = "../inputs/info.txt"
info = loadtxt(dataFile,comments="%")
n_S = int(info[0])
n_s = int(info[1])
n_phis_cal = int(info[2])
n_phis_val = int(info[3])
n_times = int(info[4])
var = float(info[5])
inad_type = int(info[6])

n_phis_tot = n_phis_cal + n_phis_val
dataFile ='gammas';

d = loadtxt(dataFile,comments="%")

filename = 'gammas-squared/'+str(n_S)+'-'+str(n_s)+'-'+str(n_phis_cal)+'-'+str(n_phis_val)+'-'+str(n_times)+'-'+str(var)+'-'+str(inad_type);

file = open(filename,'a')

d_cal = d[:(n_s*n_times*n_phis_cal)]
d_val = d[(n_s*n_times*n_phis_cal):]
# NEEDS TO SEPARATE CAL AND VAL
n_d_cal = n_s * n_times * n_phis_cal
n_d_val = n_s * n_times * n_phis_val
n_g5c = (d_cal < .05).sum()
n_g1c = (d_cal < .01).sum()
frac_g5c = float(n_g5c)/n_d_cal
frac_g1c = float(n_g1c)/n_d_cal
bon_g5c = (d_cal < (.05/n_d_cal)).sum()
bon_g1c = (d_cal < (.01/n_d_cal)).sum()

n_g5v = (d_val < .05).sum()
n_g1v = (d_val < .01).sum()
frac_g5v = float(n_g5v)/n_d_val
frac_g1v = float(n_g1v)/n_d_val
bon_g5v = (d_val < (.05/n_d_val)).sum()
bon_g1v = (d_val < (.01/n_d_val)).sum()
file.write(str(frac_g5c)+' '+str(bon_g5c)+' '+str(frac_g1c)+' '+str(bon_g1c)+' '+str(frac_g5v)+' '+str(bon_g5v)+' '+str(frac_g1v)+' '+str(bon_g1v)+'\n')
file.close()
# if gamma < .004:
#     file.write(str(gamma))
# else:
#     file.write(str(round(gamma,2)))
# file.write('\n')
