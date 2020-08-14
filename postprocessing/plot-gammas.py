import matplotlib
from matplotlib import rc
from matplotlib import rcParams
import matplotlib.pyplot as plt
import numpy as np
from numpy import loadtxt
from scipy.interpolate import interp1d
from operator import add
from operator import sub

rc('text',usetex=True)
rcParams.update({'figure.autolayout': True})
font={'family' : 'Times',
    'weight' : 'normal',
    'size' :18}
matplotlib.rc('font',**font)

max_s = 4
S = 20
averages = np.zeros((max_s,8))

for i in range(1,max_s+1):
    #dataFile = "gammas-collected/"+str(S)+"-"+str(i)+"-3-3-10-0.001-1"
    dataFile = "gammas-squared/"+str(S)+"-"+str(i)+"-3-3-10-0.001-2"
    d = loadtxt(dataFile,comments="%")
    ad = np.mean(d,axis=0)
    # d[:,0] = frac5c
    # d[:,1] = bon5c
    # d[:,2] = frac1c
    # d[:,3] = bon1c
    # d[:,4] = frac5v
    # d[:,5] = bon5v
    # d[:,6] = frac1v
    # d[:,7] = bon1v
    averages[(i-1),:] = ad

plt.figure(1)
x = range(1,max_s+1)
plt.plot(x, averages[:,0],linewidth=3,label='Calibration')
plt.plot(x, averages[:,4],linewidth=3,label='Validation')
plt.axhline(y=0.05,linewidth=1.5,color='black',alpha=.3,linestyle='--')
plt.xlabel('$s$');# (number of species included in calibration)')
#plt.ylabel('Average fraction of $\gamma$-values below 0.05')
plt.ylabel('$f_\gamma( 10,s, p, 100, 0.05)$')
plt.legend(loc=0)
#plt.savefig('/Users/rebeccam/Documents/papers/enriched-glv/rawfigs/S''%s''gamma05.pdf' %(S))

plt.figure(2)
x = range(1,max_s+1)
plt.plot(x, averages[:,2],linewidth=3,label='Calibration')
plt.plot(x, averages[:,6],linewidth=3,label='Validation')
plt.axhline(y=0.01,linewidth=1.5,color='black',alpha=.3,linestyle='--')
plt.xlabel('$s$');# (number of species included in calibration)')
#plt.ylabel('Average fraction of $\gamma$-values below 0.01')
plt.ylabel('$f_\gamma( 10,s, p, 100, 0.01)$')
plt.legend(loc=0)
#plt.savefig('/Users/rebeccam/Documents/papers/enriched-glv/rawfigs/S''%s''gamma01.pdf' %(S))
plt.show()

    
# dataFile = "/point/info.txt"
# info = loadtxt(dataFile,comments="%")
# n_S = int(info[0])
# n_s = int(info[1])
# n_phis_cal = int(info[2])
# n_phis_val = int(info[3])
# n_times = int(info[4])
# var = int(info[5])
# inad_type = int(info[6])

  
    # plt.figure(k*n_s+i)
    # plt.fill_between(times,r2,r3,facecolor='blue',alpha=.3)
    # plt.fill_between(times,r1,r4,facecolor='blue',alpha=.1)
    # plt.plot(times,rmean,'b',linewidth=2,label='Augmented model')
    # plt.plot(times,state,'r^',markersize=7,label='Detailed model (data)')
    # # print(datatimes)
    # # print(times)
    # red = interp1d(datatimes,state_red)
    # plt.plot(times,red(times),'g', linewidth=2,label='Reduced model')
  
    # plt.xlabel('Time')
    # # plt.ylabel(spec_names[i]+ylabel )
    # if k < n_phis_cal:
    #   plt.title('Calibration scenario '+str(k+1)+' of '+str(n_phis_cal))
    # else:
    #   plt.title('Prediction scenario '+str(k+1-n_phis_cal)+' of '+str(n_phis_val))

    # plt.ylabel('$x_{}$'.format(i+1),fontsize=28)
    # plt.xlim(times[0],times[-1])
    # plt.locator_params(nbins=5)
    # plt.legend(loc=0)
    # # plt.show()
   # # plt.savefig('../../../../documents/juq15/rawfigs/h2-op/smooth-'
   # #     '%s' '-' '%s''.pdf' %(k,i))
    # # plt.savefig('../../talks/2018/ini/figures/smooth-'
    # #    'S' '%s' '-s' '%s' '-phi' '%s' '-spec' '%s''.pdf' %(n_S,n_s,k,i))
    # # plt.savefig('red-plots/smooth-'
    # #         '%s' '-' '%s''.pdf' %(k,i))
    # # plt.savefig('red-plots/smooth-pred-'
    # #        '%s' '-' '%s''.pdf' %(k,i))

# plt.show()
