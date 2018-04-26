import numpy as np 
import pandas as pd

import matplotlib.pyplot as plt 
import seaborn as sns 

sns.set()

obd_df =  np.loadtxt("../data/fuckyou.csv", delimiter=",")


plt.hist(x = obd_df[:,1], bins = obd_df[:,0])
plt.xlabel(r"radius $r$")
plt.ylabel(r"density $\rho(r)$")  
plt.show()


