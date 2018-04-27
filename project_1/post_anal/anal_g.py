import numpy as np 
import pandas as pd

import matplotlib.pyplot as plt 
import seaborn as sns 

sns.set()


ni_filename = "../data/obd_NM_NIA_a_0.400000_b_1.000000_step_0.100000_np_100_nd_3.csv"  
i_filename = "../data/obd_NM_IA_a_0.400000_b_1.000000_step_0.100000_np_100_nd_3.csv" 

ni_obd_df =  np.loadtxt(ni_filename, delimiter=",",  skiprows = 1)
i_obd_df =  np.loadtxt(i_filename, delimiter=",",  skiprows = 1)


plt.hist(ni_obd_df[:,1], bins = ni_obd_df[:,0], label = "non interactive", alpha = 0.6)
plt.hist(i_obd_df[:,1], bins = i_obd_df[:,0], label = "interactive", alpha = 0.6)

plt.xlabel(r" $r$")
plt.ylabel(r"$\rho(r)$")  

plt.legend()
plt.show()


