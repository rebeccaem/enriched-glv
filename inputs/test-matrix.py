import matplotlib
import numpy as np
from numpy import loadtxt

n_S = 2
n_s = 1
n_phis_cal = 1
n_phis_val = 0
n_times = 20
var = .001
inad_type = 1
filename = 'inputs/info.txt';
file = open(filename,'w')
file.write(str(n_S)+'\n'+str(n_s)+'\n'+str(n_phis_cal)+'\n'+str(n_phis_val)+'\n'+str(n_times)+'\n'+str(var)+'\n'+str(inad_type)+'\n')
file.close()

filename = 'inputs/matrix.txt';
file = open(filename,'w')

# b = np.random.lognormal(0,1,n_S)
# B = np.random.lognormal(0,1,(n_S,n_S))
# B = np.tril(B,-1)
# B = B + B.T
# C = np.sqrt(n_S)*np.diag(np.random.lognormal(0,1,n_S))
# # C = n_S*np.diag(np.random.lognormal(0,1,n_S))
# A = -(B + C)
#A = np.array([[-3,-1],[-1,-2]])
A = np.array([[-3,-1],[-1,-2]])
b = np.array([5,3])
for i in range(n_S):
  file.write(str(b[i])+'\n')
for i in range(n_S):
  for j in range(n_S):
    file.write(str(A[i][j])+'\n')

file.close()
