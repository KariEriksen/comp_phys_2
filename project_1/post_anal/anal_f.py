import numpy as np 
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt 
import seaborn as sns

sns.set()


f1 = "../data/GD_NM_IA_np_5_nd_3_sa_0.400000.csv"
f2 = "../data/GD_NM_IA_np_5_nd_3_sa_0.300000.csv"
f3 = "../data/GD_NM_IA_np_5_nd_3_sa_0.600000.csv"

fn = [f1, f2, f3]
arrays = []

for fname in fn:
	arrays.append(np.loadtxt(fname, delimiter = ","))
