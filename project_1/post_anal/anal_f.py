import numpy as np 
import pandas as pd
import numpy as np

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt 
import seaborn as sns

sns.set()


f1 = "../data/GD_NM_IA_np_5_nd_3_sa_0.400000.csv"
f2 = "../data/GD_NM_IA_np_5_nd_3_sa_0.450000.csv"
f3 = "../data/GD_NM_IA_np_5_nd_3_sa_0.600000.csv"

fn = [f1, f2, f3]
alphas = [0.4, 0.45, 0.6]
arrays = []

for fname in fn:
	arrays.append(np.loadtxt(fname, delimiter = ","))


fig = plt.figure()
ax = fig.gca(projection='3d') 		

for i in range(len(arrays)):
	array = arrays[i]
	alpha_start = alphas[i]
	
	alpha  = array[:,0]
	el = array[:,1]
	iterat = np.arange(len(el))

	alpha_str = "{:.2f}".format(alpha_start)
	plt.plot(alpha, el, iterat, label = r"start $\alpha = $ " + alpha_str, alpha = 0.8) 


minima_label = r"$\alpha_{min} = $"+"{:.4e}".format(alpha[-1])
ax.text(alpha[-1], el[-1],  10, minima_label, fontsize = "smaller")

plt.title(r"Gradient descent optimization of variational parameter $\alpha$", y = 1.08)

ax.set_xlabel(r"$\alpha$")
ax.set_zlabel("i ")
ax.set_ylabel(r"$\langle E_L (\alpha ) \rangle$")

for tick in ax.get_yticklabels():
    tick.set_rotation(90)

ax.yaxis.labelpad = 20

plt.legend(loc = "best")
plt.show()