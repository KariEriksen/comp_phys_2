"""
Analysis of data produced for project 2 task c.
Does exactly the same as anal_b.py, but extends the functionality
to provide comparison of importance sampling and naive metropolis algorithm.
Comparison needs:
    Runtimes for naive_mh and importance
        - should be available from NM_NQS and IM_NQS meta files
    Variance in E_l for both variants
        -   should be easily obtainable from already imported data
"""


import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
import re
import seaborn as sns
from numpy import log2, zeros, mean, var, sum, loadtxt, arange, array, cumsum, dot, transpose, diagonal, floor
from numpy.linalg import inv
from time import time
import numpy as np

sns.set()
sns.set_context("paper")

def block(x):
    # preliminaries
    n = len(x); d = int(log2(n));
    s, gamma = zeros(d), zeros(d);
    mu = mean(x); t0 = time()

    # estimate the auto-covariance and variances 
    # for each blocking transformation
    for i in arange(0,d):
        n = len(x)
        # estimate autocovariance of x
        gamma[i] = (n)**(-1)*sum( (x[0:(n-1)]-mu)*(x[1:n]-mu) )
        # estimate variance of x
        s[i] = var(x)
        # perform blocking transformation
        x = 0.5*(x[0::2] + x[1::2])
   
    # generate the test observator M_k from the theorem
    M = (cumsum( ((gamma/s)**2*2**arange(1,d+1)[::-1])[::-1] )  )[::-1]

    # we need a list of magic numbers
    q =array([6.634897,9.210340, 11.344867, 13.276704, 15.086272, 16.811894,
        18.475307, 20.090235, 21.665994, 23.209251, 24.724970, 26.216967,
        27.688250, 29.141238, 30.577914, 31.999927, 33.408664, 34.805306,
        36.190869, 37.566235, 38.932173, 40.289360, 41.638398, 42.979820,
        44.314105, 45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

    # use magic to determine when we should have stopped blocking
    for k in arange(0,d):
        if(M[k] < q[k]):
            break
    if (k >= d-1):
        print("Warning: Use more data")
    ans = s[k]/2**(d-k)
    #print( "Runtime: %g sec" % (time()-t0)); print( "Blocking Statistics :")
    #print("average            iterations      std. error")
    #print ("%8g %20g %15g" % (mu, k, ans**.5))
    return ans

filenames = sorted(np.array(os.listdir("../data/b_data")))
filenames = [name for name in filenames if name != "dummy" and name != "time_iter.csv"]

# Perform blocking on results
gamma_vals = [0.01, 0.26, 0.51, 0.76]
blocking_data = []
energies = []
for gamma in gamma_vals:
    tmp_block = []
    tmp_energy = []
    for filename in filenames: 
        if str(gamma) in filename:
            t_filename = "../data/b_data/"+filename
            A = loadtxt(t_filename)
            tmp_block.append(block(A))
            tmp_energy.append(np.mean(A))

    blocking_data.append(tmp_block)
    energies.append(tmp_energy)


time_data = np.loadtxt("../data/b_data/time_iter.csv")
mean_time = np.mean(np.array(time_data))

# Plotting
x = range(len(energies[0]))
fig, ax = plt.subplots(figsize=(9, 7))
for i in range(len(energies)):
    plt.errorbar(x, energies[i], fmt = "^-",
            barsabove = True,
            yerr = np.array(blocking_data[i], dtype = float),
            label = r"$\gamma = ${:2f}".format(gamma_vals[i]))

plt.xticks(x, rotation = 45, size = "medium")
plt.legend()
plt.xlabel(r"Iteration")
plt.ylabel(r"$\langle E \rangle $")
plt.title("Importance sampling with varied gamma | mean time = {:.2g}s".format(mean_time))

#plt.savefig("../report/figures/importance")
plt.show()
#plt.clf()

