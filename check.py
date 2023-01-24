import numpy as np
import matplotlib.pyplot as plt
import inputs
import amp as amp

model = inputs.model
ch = inputs.ch

data1 = np.genfromtxt(inputs.cnt_file+"_fp.dat").T
data2 = np.genfromtxt(inputs.cnt_file+"_sp.dat").T

plt.plot(data1[0], data1[1])
plt.scatter(data1[0], data1[2])
plt.scatter(data1[0], data1[1])
plt.plot(data1[0], data1[2])

plt.plot(data2[0], data2[1])
plt.plot(data2[0], data2[2])

plt.yscale("log")
plt.xscale("log")

plt.show()
plt.close()

