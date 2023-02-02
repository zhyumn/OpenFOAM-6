import numpy as np
import os
import matplotlib.pyplot as plt
folder = os.path.split(os.path.realpath(__file__))[0]
file = open(folder+ "/postProcessing/cellMax(mag(subtract(T_ex,T)))/0/volFieldValue.dat", 'r')
data = [s.split("\t") for s in file.read().split("\n")[4:-1]]
data = np.array([[float(s2) for s2 in s1] for s1 in data])

file = open(folder+ "/postProcessing/volAverage(mag(subtract(T_ex,T)))/0/volFieldValue.dat", 'r')
data2 = [s.split("\t") for s in file.read().split("\n")[4:-1]]
data2 = np.array([[float(s2) for s2 in s1] for s1 in data2])



tT=3e0 
line_T_max, = plt.plot(data[:,0], data[:,1]/tT, "--",color= "b")
line_T_mean, =plt.plot(data2[:,0], data2[:,1]/tT,color= "b")
plt.yscale('log')


tp=3e5

folder = os.path.split(os.path.realpath(__file__))[0]
file = open(folder+ "/postProcessing/cellMax(mag(subtract(p_ex,p)))/0/volFieldValue.dat", 'r')
data = [s.split("\t") for s in file.read().split("\n")[4:-1]]
data = np.array([[float(s2) for s2 in s1] for s1 in data])

file = open(folder+ "/postProcessing/volAverage(mag(subtract(p_ex,p)))/0/volFieldValue.dat", 'r')
data2 = [s.split("\t") for s in file.read().split("\n")[4:-1]]
data2 = np.array([[float(s2) for s2 in s1] for s1 in data2])

line_p_max, =plt.plot(data[:,0], data[:,1]/tp, "--", color= "r")
line_p_mean, =plt.plot(data2[:,0], data2[:,1]/tp,color= "r")
plt.legend([line_T_mean, line_p_mean], ['Temperature', 'Pressure'],loc='lower right',frameon=False,fontsize=13)
plt.xlabel("time/s", fontsize=12)
plt.axhline(y = 1, color='black', linestyle='dotted')
plt.text(0.2e-6,0.2,"Maximum error", fontsize=14)
plt.text(0.8e-6,2e-4,"Average error", fontsize=14)
plt.axhline(y = 1, color = 'black', linestyle = 'dotted')
print("fig saved")
plt.savefig(folder + "//error.png", dpi=400 , bbox_inches='tight')