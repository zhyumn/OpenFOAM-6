import numpy as np
import os
import matplotlib.pyplot as plt
folder = os.path.split(os.path.realpath(__file__))[0]
file = open(folder+ "/postProcessing/cellMax(mag(subtract(T_ex,T)))/2e-06/volFieldValue.dat", 'r')
data = [s.split("\t") for s in file.read().split("\n")[4:-1]]
data = np.array([[float(s2) for s2 in s1] for s1 in data])

file = open(folder+ "/postProcessing/volAverage(mag(subtract(T_ex,T)))/2e-06/volFieldValue.dat", 'r')
data2 = [s.split("\t") for s in file.read().split("\n")[4:-1]]
data2 = np.array([[float(s2) for s2 in s1] for s1 in data2])



tT=1e0 
plt.plot(data[:,0], data[:,1]/tT)
plt.plot(data2[:,0], data2[:,1]/tT)
plt.yscale('log')


tp=1e5

folder = os.path.split(os.path.realpath(__file__))[0]
file = open(folder+ "/postProcessing/cellMax(mag(subtract(p_ex,p)))/2e-06/volFieldValue.dat", 'r')
data = [s.split("\t") for s in file.read().split("\n")[4:-1]]
data = np.array([[float(s2) for s2 in s1] for s1 in data])

file = open(folder+ "/postProcessing/volAverage(mag(subtract(p_ex,p)))/2e-06/volFieldValue.dat", 'r')
data2 = [s.split("\t") for s in file.read().split("\n")[4:-1]]
data2 = np.array([[float(s2) for s2 in s1] for s1 in data2])

plt.plot(data[:,0], data[:,1]/tp)
plt.plot(data2[:,0], data2[:,1]/tp)

plt.axhline(y = 1, color = 'r', linestyle = '--')
print("fig saved")
plt.savefig(folder + "//error.png", dpi=400 , bbox_inches='tight')