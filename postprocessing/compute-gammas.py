import matplotlib
import numpy as np
from numpy import loadtxt
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

burnin = 200;
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
dataFile = "../inputs/datafile.txt"
d = loadtxt(dataFile,comments="%")
dataFile = "sfp_qoi_seq.dat"
q = loadtxt(dataFile,comments="%")
q = q[burnin:,:]
# q = q.flatten()

sigma = np.sqrt(var)
# print(sigma)

filename = 'gammas';
file = open(filename,'w')

# for i in range(n_phis_tot*n_times*n_s):
for i in range(n_phis_tot):
    for j in range(n_times):
        for k in range(n_s):
            d_idx = i*n_times*(n_S) + j*(n_S) + k
            # print("idx = ",d_idx)
            q_idx = i*n_times*2*(n_s) + j*2*(n_s) + k
            # (a, b) = divmod(i,n_s)
            # j = a*n_s + b
            f = [];
            # for j in range(len(q[:,0])):
            for l in range(q.shape[0]):
                mu = q[l,q_idx];
                f.extend(np.random.normal(mu,sigma,100));

            density = gaussian_kde(f)
            # print(d[d_idx,2])
            pofy = density.evaluate(d[d_idx,2])
            minf = np.min(f)
            maxf = np.max(f)
            xs = np.linspace(minf,maxf,200)
            fpdf = density.evaluate(xs)
            # plt.figure()
            # plt.plot(xs,fpdf)
            # plt.show()
            for l in range(len(fpdf)):
                #print h[j]
                if fpdf[l] > pofy:
                    fpdf[l] = 0.
            gamma = np.trapz(fpdf,xs);
            # print(gamma)
            if gamma < .004:
                file.write(str(gamma))
            else:
                file.write(str(round(gamma,2)))
            file.write('\n')
file.close()
