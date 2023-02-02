import numpy as np
import os
import matplotlib.pyplot as plt
folder = os.path.split(os.path.realpath(__file__))[0]
file = open(folder+ "/ISAT_VLE/VLEtime", 'r')
data = [s.split(",") for s in file.read().split("\n")[:-1]]
data = np.array([[float(s2) for s2 in s1] for s1 in data])

file = open(folder+ "/ISAT_VLE/cputotal", 'r')
data2 = [s.split(",") for s in file.read().split("\n")[:-1]]
data2 = np.array([[float(s2) for s2 in s1] for s1 in data2])




def moving_average(interval,windowsize):
    window = np.ones(int(windowsize))/float(windowsize)
    re = np.convolve(interval,window,'same')
    return re
#----------------------------------------------


plt.plot(data[:,0][0::10], moving_average(data[:,1],10)[0::10],  color= "r")
plt.plot(data2[:,0][0::10], moving_average(data2[:,1]-data[:,1],10)[0::10],'--', color= "b")

plt.xlim(0,2e-6)
plt.ylim(0,3.5)
plt.legend( ['VLE time', 'Others'],loc='upper left',frameon=False,fontsize=13)
plt.xlabel("physical time (s)", fontsize=12)
plt.ylabel("cpu time (s)", fontsize=12)

plt.savefig(folder + "//time_noISAT.png", dpi=400 , bbox_inches='tight')

#----------------------------------------------


