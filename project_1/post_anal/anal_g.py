import numpy as np 
import pandas as pd

import matplotlib.pyplot as plt 
import seaborn as sns 

sns.set()

obd_df =  np.loadtxt("../data/obd_test.csv", delimiter=",",  skiprows = 1)



plt.hist(obd_df[:,1], bins = obd_df[:,0])
plt.xlabel(r" $r$")
plt.ylabel(r"$\rho(r)$")  
plt.show()


