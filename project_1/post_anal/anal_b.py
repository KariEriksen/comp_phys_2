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
    q =array([6.634897,9.210340, 11.344867, 13.276704, 15.086272, 16.811894, 18.475307, 20.090235, 21.665994, 23.209251, 24.724970, 26.216967, 27.688250, 29.141238, 30.577914, 31.999927, 33.408664, 34.805306, 36.190869, 37.566235, 38.932173, 40.289360, 41.638398, 42.979820, 44.314105, 45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

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

filenames = np.array(os.listdir("../data/temp_data"))

#With error plot
alpha_regex = r".+_a_(0.\d{6}).+"
alpha_regex_obj = re.compile(alpha_regex)

header_regex =r"# N_p:(\s\d+)\| N_d:(\s\d+)\| N_mc:(\s\d+)"
header_regex_obj = re.compile(header_regex)

mask = [not f.endswith("data.csv") for f in filenames]
filenames = sorted(filenames[mask], key = lambda x: float(alpha_regex_obj.search(x).group(1)))
stds = [0, 0]

N_ps = (1, 10, 100, 500)
N_mc = 0
N_ds = (1, 2, 3)
sim_types = [v for v in sys.argv[1:]]

fast = False

for N_d_i in N_ds:
    for N_p_i in N_ps:
        mean_n_times = pd.read_csv("../data/temp_data/"+
                sim_types[0]+
                "_"+sim_types[1][3:]+
                "_meta"+
                "_np_"+str(N_p_i)+
                "_nd_"+str(N_d_i)+
                "_data.csv")
        
        for i in range(len(sim_types)):
            sim_type = sim_types[i]
            tmp = []
            for filename in filenames: 
                if filename.startswith(sim_type):
                    t_filename = "../data/temp_data/"+filename
                    
                    f_o = open(t_filename)
                    header = f_o.readline()
                    f_o.close()

                    match_o = header_regex_obj.search(header)
                    N_p = match_o.group(1)
                    N_d = match_o.group(2)
                    N_mc = match_o.group(3)

                    if fast:
                        break 
                    elif int(N_p) == N_p_i and int(N_d) == N_d_i:
                        A = loadtxt(t_filename, delimiter = "\n", skiprows = 1)
                        tmp.append(block(A))
            stds[i] = tmp
            

        fig, ax = plt.subplots(figsize=(9, 7))

        time_anal = np.average(mean_n_times["analytic_time"])
        try:
            time_num = np.average(mean_n_times["numeric_time"])
            num_min_index = np.where(
                    mean_n_times.numeric_energy == np.min(mean_n_times.numeric_energy))
            num_min_cord = (mean_n_times.alpha[num_min_index[0][0]],
                mean_n_times.numeric_energy[num_min_index[0][0]])

            t_r = time_num / time_anal
        except IndexError:
            pass 


        anal_min_index = np.where(
                mean_n_times.analytic_energy == np.min(mean_n_times.analytic_energy))
        anal_min_cord = (mean_n_times.alpha[anal_min_index[0][0]],
                mean_n_times.analytic_energy[anal_min_index[0][0]])

        title =  sim_types[0]+"_"+sim_types[1][3:]+ "_" + "VMC with "
        title += r"$N_P = ${} | $N_D = ${} | $N_{{MC}}$ = {:.0e}".format(N_p_i, N_d_i, int(N_mc))
        title += r"| time ratio : $\frac{T_N}{T_A} = $" 
        title += "{:.2e}".format(t_r)


        if fast:
            mean_n_times.plot("alpha", ["analytic_energy", "numeric_energy"],
                    legend = True,
                    title = title,
                    ax = ax) 
        else:
            x = mean_n_times.alpha.values
            num = mean_n_times.numeric_energy.values
            anal = mean_n_times.analytic_energy.values

            plt.errorbar(x, num, fmt = "x-",
                    barsabove = True,
                    yerr = np.array(stds[1], dtype = float),
                    label = "numeric_energy")
            plt.errorbar(x, anal, fmt = "^-",
                    barsabove = True,
                    yerr = np.array(stds[0], dtype = float),
                    label = "analytic_energy") 
            plt.title(title)
            plt.xticks(x, rotation = 45, size = "medium")
            plt.legend()

        a_min_txt = r"Analytic GS @ $\alpha = ${} with $\langle E \rangle / N_P =$ {:.2f}".format(anal_min_cord[0], anal_min_cord[1]/float(N_p_i))

        n_min_txt = r"Numeric GS @ $\alpha = ${} with $\langle E \rangle / N_P =$ {:.2f}".format(num_min_cord[0], num_min_cord[1]/float(N_p_i))

        ax.annotate(a_min_txt, xy=anal_min_cord,  xycoords='data',
                xytext=(50, 120), textcoords='offset points',
                arrowprops=dict(arrowstyle="->",
                                connectionstyle="arc3,rad=0.2"))


        ax.annotate(n_min_txt, xy=num_min_cord,  xycoords='data',
                xytext=(-60, 150), textcoords='offset points',
                arrowprops=dict(arrowstyle="->",
                                connectionstyle="arc3,rad=0.2"))        
        plt.xlabel(r"$\alpha$")
        plt.ylabel(r"$\langle E \rangle $")

        figname = sim_types[0]+"_"+sim_types[1][3:] 
        figname += "_np_" + str(N_p_i) + "_nd_" + str(N_d_i) + ".pdf" 
        plt.savefig("../report/figures/" + figname)
        plt.clf()
