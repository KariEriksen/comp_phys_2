import pandas as pd 
import matplotlib.pyplot as plt
import seaborn

filename = "../data/vmc_N_10_D_2.txt"
df = pd.read_csv(filename, skiprows = 0, header = 1, index_col = False)

df.plot(x = " alpha", y = "mean(E_l)")
plt.show()
