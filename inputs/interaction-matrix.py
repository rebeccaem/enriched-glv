import matplotlib
import numpy as np
from numpy import loadtxt
# np.random.seed(21)

n_S = 20
n_s = 4
n_phis_cal = 3
n_phis_val = 3
n_times = 10
var = .001
#model inadequacy type: 1 = absolute value, 2 = squared, 3 = nonlinear diagonal, 4 = nonlinear mixed
#10 = memory, 11 = algebraic
inad_type = 1
filename = 'inputs/info.txt';
file = open(filename,'w')
file.write(str(n_S)+'\n'+str(n_s)+'\n'+str(n_phis_cal)+'\n'+str(n_phis_val)+'\n'+str(n_times)+'\n'+str(var)+'\n'+str(inad_type)+'\n')
file.close()

# choose if A is diagonally dominant or not:
# set to 0 for not, set to 1 to be dd
c = 1;

filename = 'inputs/matrix.txt';
file = open(filename,'w')

# b = np.sqrt(n_S)*np.random.lognormal(0,1,n_S)
# b = n_S*[n_S]
B = np.random.lognormal(0,1,(n_S,n_S))
B = np.tril(B,-1)
B = B + B.T
offdiag = np.sum(B,axis=1)
# C = np.sqrt(n_S)*np.diag(np.random.lognormal(0,1,n_S))
if c == 0:
    C = np.diag(np.random.lognormal(0,n_S,n_S))
elif c == 1:
    # make A diagonally dominant with below
    C = np.diag(offdiag + np.random.lognormal(0,1,n_S))

print(np.max(C))
A = -(B + C)
b = np.ones((n_S,1))
b = np.max(C)*b

for i in range(n_S):
  file.write(str(b[i][0])+'\n')
for i in range(n_S):
  for j in range(n_S):
    file.write(str(A[i][j])+'\n')

file.close()
