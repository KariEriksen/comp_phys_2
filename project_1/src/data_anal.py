import matplotlib.pyplot as plt
import os
import re
import seaborn
from numpy import log2, zeros, mean, var, sum, loadtxt, arange, array, cumsum, dot, transpose, diagonal, floor
from numpy.linalg import inv
from time import time
import numpy as np

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
    q =array([6.634897,9.210340, 11.344867, 13.276704, 15.086272, 16.811894, 18.475307, 20.090235, 21.665994, 23.209251, 24.724970, 26.216967, 27.688250, 29.141238, 30.577914, 31.999927, 33.408664, 34.805306, 36.190869, 37.566235, 38.932173, 40.289360, 41.638398, 42.979820, 44.314105, 45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

    # use magic to determine when we should have stopped blocking
    for k in arange(0,d):
        if(M[k] < q[k]):
            break
    if (k >= d-1):
        print("Warning: Use more data")
    ans = s[k]/2**(d-k)
    print( "Runtime: %g sec" % (time()-t0)); print( "Blocking Statistics :")
    print("average            iterations      std. error")
    print ("%8g %20g %15g" % (mu, k, ans**.5))
    return ans

filenames = os.listdir("../data/")

alpha_regex = r"\w+_a_(0.\d{6}).+"
alpha_regex_obj = re.compile(alpha_regex)

header_regex =r"# N_p:(\s\d+)\| N_d:(\s\d+)\| N_mc:(\s\d+)"
header_regex_obj = re.compile(header_regex)

filenames = sorted(filenames, key = lambda x: float(alpha_regex_obj.search(x).group(1)))
means = []
stds = []
times = []
x = np.arange(0.3, 0.8, 0.05)

N_p = 0
N_mc = 0
N_d = 0
sim_type = "NM_NIA"

for filename in filenames: 
    if filename.startswith(sim_type):
        t_filename = "../data/"+filename
        
        f_o = open(t_filename)
        header = f_o.readline()
        time_mc = f_o.readline()
        f_o.close()

        times.append(float(time_mc))
        if N_p == 0:
            match_o = header_regex_obj.search(header)
            N_p = match_o.group(1)
            N_d = match_o.group(2)
            N_mc = match_o.group(3)

        A = loadtxt(t_filename, delimiter = "\n", skiprows = 2)
        stds.append(block(A))
        means.append(mean(A))
    else:
        print("fuck you")

means = np.array(means)/float(N_p)
plt.errorbar(x, means, yerr = stds)
plt.title(sim_type + "_"+ r"VMC with $N_P = ${} | $N_D = ${}| average time[s]: {:g}".format(N_p, N_d, 
    np.average(np.array(times))))
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$E_L∕N_p$")
plt.show()
