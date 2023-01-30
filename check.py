import numpy as np
import matplotlib.pyplot as plt
import inputs

model = inputs.model
ch = inputs.ch

data1 = np.genfromtxt(inputs.cnt_file+"_fp.dat").T
data2u = np.genfromtxt(inputs.cnt_file+"_sp_u.dat").T
data2d = np.genfromtxt(inputs.cnt_file+"_sp_d.dat").T
data3u = np.genfromtxt(inputs.cnt_file+"_tp_u.dat").T
data3d = np.genfromtxt(inputs.cnt_file+"_tp_d.dat").T
data4 = np.genfromtxt(inputs.cnt_file+"_op.dat").T

plt.plot(data1[0], data1[1], c='red')
plt.plot(data1[0], data1[2], c='red')

plt.plot(data2u[0], data2u[1], c='blue')
plt.plot(data2d[0], data2d[1], c='blue')

plt.plot(data3u[0], data3u[1], c='orange')
plt.plot(data3d[0], data3d[1], c='orange')
#   plt.plot(data3[0], data3[1], c='green')
#   plt.plot(data3[0], data3[2], c='green')

plt.yscale("log")
plt.xscale("log")

plt.show()
plt.close()
