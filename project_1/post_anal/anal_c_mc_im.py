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



def plotter_func(A, N_d, N_p, N_mc, alphas,sim_type):

    title = sim_type + " " + str(int(N_d)) + "D | " + str(int(N_p)) + r"P | $\alpha$ " + str(alphas)
    cycles = np.arange(int(N_mc))

    plt.plot(cycles, A)
    figname = sim_type + "_" + str(int(N_d)) + "D_" + str(int(N_p)) + "P"
    
    plt.xlabel("MC cycles")
    plt.ylabel(r"$E_L$")
    plt.title(title)
    plt.savefig(figname+".pdf")
    plt.clf()
    plt.cla()
    return


filenames = np.array(os.listdir("../data/"))

#With error plot
alpha_regex = r".+_a_(0.\d{6}).+"
alpha_regex_obj = re.compile(alpha_regex)

header_regex =r"# N_p:(\s\d+)\| N_d:(\s\d+)\| N_mc:(\s\d+)"
header_regex_obj = re.compile(header_regex)

mask = [not f.endswith("data.csv") for f in filenames]
filenames = sorted(filenames[mask], key = lambda x: float(alpha_regex_obj.search(x).group(1)))
stds = [0, 0]

N_ps = (500,)
N_mc = 0
N_ds = (3,)
sim_types = [v for v in sys.argv[1:]]

for N_d_i in N_ds:
    for N_p_i in N_ps:
        mean_n_times = pd.read_csv("../data/IM_NIA_meta_np_500_nd_3_dt_0.100000_data.csv")
        
        for i in range(len(sim_types)):
            sim_type = sim_types[i]
            tmp = []
            t_filename = "../data/IM_NIA_a_0.400000_b_1.000000_step_1.000000_np_500_nd_3_dt_0.100000.csv"
            
            f_o = open(t_filename)
            header = f_o.readline()
            f_o.close()

            match_o = header_regex_obj.search(header)
            N_p = match_o.group(1)
            N_d = match_o.group(2)
            N_mc = match_o.group(3)

            A = loadtxt(t_filename, delimiter = "\n", skiprows = 1)
            alpha = float(re.match(alpha_regex,t_filename).group(1))
            plotter_func(A, N_d, N_p, N_mc, alpha,sim_type)

