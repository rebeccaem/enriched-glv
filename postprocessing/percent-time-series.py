import matplotlib
from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
from numpy import loadtxt
from scipy.interpolate import interp1d
from operator import add
from operator import sub

rc('text',usetex=True)
font={'family' : 'normal',
    'weight' : 'normal',
    'size' :14}
matplotlib.rc('font',**font)

dataFile = "../inputs/info.txt"
info = loadtxt(dataFile,comments="%")
n_S = int(info[0])
n_s = int(info[1])
n_phis_cal = int(info[2])
n_phis_val = int(info[3])
n_times = int(info[4])
var = int(info[5])
inad_type = int(info[6])

n_phis_tot = n_phis_cal + n_phis_val
#dataFile = "datafile-det-pred.txt"
dataFile = "../inputs/datafile.txt"
# dataFile = "datafile-det-pred-dense.txt"
data = loadtxt(dataFile,comments="%")
datatimes = data[0:(n_times*n_S):n_S,1]; #all times are the same
times = np.linspace(0,datatimes[-1],num=n_times,endpoint=True)

dataFile = "../inputs/datafile-reduced.txt"
data_red = loadtxt(dataFile,comments="%")

dataFile = "qoi-stats"
q = loadtxt(dataFile,comments="%")

# spec_names = ('H$_2$', 'O$_2$', 'H', 'O', 'OH', 'HO$_2$', 'H$_2$O', 'H$_2$O$_2$',
#     'N$_2$', 'Temperature')

sigma = np.sqrt(var)
for k in range(0,n_phis_tot):
#for k in range(0,1):
#for k in range(4,5):
  for i in range(0,n_s): 
    state = data[(k*n_times*n_S + i):((k+1)*n_times*n_S + i):n_S,2];
    state_red = data_red[(k*n_times*n_s + i):((k+1)*n_times*n_s + i):n_s,2];
    # err = np.random.normal(0,sigma,n_times)
    # state = state + err
    # state[np.where(state<0)] = 0
    # ylabel = ' [mol/m$^3$]'
    rmean = [];
    r1 = []; r2 = []; r3 = []; r4 = [];
    for j in range(0,n_times):
      rmean.append(q[k*2*n_s*n_times + j*2*n_s + i,0])
      r1.append(q[k*2*n_s*n_times + j*2*n_s + i,1])
      r2.append(q[k*2*n_s*n_times + j*2*n_s + i,2])
      r3.append(q[k*2*n_s*n_times + j*2*n_s + i,3])
      r4.append(q[k*2*n_s*n_times + j*2*n_s + i,4])
  
    plt.figure(k*n_s+i)
    plt.fill_between(times,r2,r3,facecolor='C0',alpha=.3)
    plt.fill_between(times,r1,r4,facecolor='C0',alpha=.1)
    #plt.plot(times,state,'-.',color='C3',linewidth=2,label='Detailed')
    plt.plot(times,state,'^',markersize=7,color='C3',label='Detailed (data)')
    # print(datatimes)
    # print(times)
    red = interp1d(datatimes,state_red)
    plt.plot(times,red(times),'--', linewidth=2,color='C1',label='Reduced')
    plt.plot(times,rmean,linewidth=2,color='C0',label='Enriched')
  
    plt.xlabel('Time')
    # plt.ylabel(spec_names[i]+ylabel )
    if k < n_phis_cal:
      # plt.title('Calibration scenario '+str(k+1)+' of '+str(n_phis_cal))
      plt.title('')
    else:
      plt.title('Validation scenario '+str(k+1-n_phis_cal)+' of '+str(n_phis_val))

    plt.ylabel('$x_{}$'.format(i+1),fontsize=28)
    plt.xlim(times[0],times[-1])
    plt.locator_params(nbins=5)
    plt.legend(loc=0)
    #plt.show()
    plt.savefig('/Users/rem/repos/documents/papers/enriched-glv/rawfigs/legend-S' '%s' '-s' '%s' '-phi' '%s' '-spec' '%s' '.pdf' %(n_S,n_s,k,i))
    #plt.savefig('/users/rebeccam/repos/documents/proposals/2019/iii-small/figs/smooth-S' '%s' '-s' '%s' '-phi' '%s' '-spec' '%s' '.pdf' %(n_S,n_s,k,i))
    # plt.savefig('red-plots/smooth-'
    #         '%s' '-' '%s''.pdf' %(k,i))
    # plt.savefig('red-plots/smooth-pred-'
    #        '%s' '-' '%s''.pdf' %(k,i))

#plt.show()
