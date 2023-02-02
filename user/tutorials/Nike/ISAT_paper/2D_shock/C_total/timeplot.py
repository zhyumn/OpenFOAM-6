import numpy as np
import os
import matplotlib.pyplot as plt
folder = os.path.split(os.path.realpath(__file__))[0]

def read_data(path):
    file = open(path+ "/VLEtime", 'r')
    data = [s.split(",") for s in file.read().split("\n")[:-1]]
    data = np.array([[float(s2) for s2 in s1] for s1 in data])
    
    file = open(path+ "/cputotal", 'r')
    data2 = [s.split(",") for s in file.read().split("\n")[:-1]]
    data2 = np.array([[float(s2) for s2 in s1] for s1 in data2])
    
    return np.array([data[:,0],data[:,1],data2[:,1]])

data_noISAT =read_data(folder+ "/time_noISAT")
data_3e_2 =read_data(folder+ "/time_3e_2")
data_1e_2 =read_data(folder+ "/time_1e_2")
data_3e_3 =read_data(folder+ "/time_3e_3")
data_1e_3 =read_data(folder+ "/time_1e_3")


    
'''
file = open(folder+ "/time_noISAT/VLEtime", 'r')
data = [s.split(",") for s in file.read().split("\n")[:-1]]
data = np.array([[float(s2) for s2 in s1] for s1 in data])

file = open(folder+ "/time_noISAT/cputotal", 'r')
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
'''
#----------------------------------------------
def moving_average(interval,windowsize):
    window = np.ones(int(windowsize))/float(windowsize)
    re = np.convolve(interval,window,'same')
    return re
#----------------------------------------------


#plt.plot(data_noISAT[0][0::10], moving_average(data_noISAT[1],10)[0::10],  color= "r")
#plt.plot(data_noISAT[0][0::10], moving_average(data_noISAT[2]-data_noISAT[1],10)[0::10],'--', color= "b")

plt.plot(data_noISAT[0][0::10], moving_average(data_noISAT[1]*0.2/(data_noISAT[2]-data_noISAT[1]),10)[0::10],  color= "r")
plt.plot(data_noISAT[0][0::10], moving_average((data_noISAT[2]-data_noISAT[1])*0.2/(data_noISAT[2]-data_noISAT[1]),10)[0::10],'--', color= "b")


plt.plot(data_1e_3[0][0::10], moving_average(data_1e_3[1]*0.2/(data_1e_3[2]-data_1e_3[1]),10)[0::10])

plt.plot(data_3e_3[0][0::10], moving_average(data_3e_3[1]*0.2/(data_3e_3[2]-data_3e_3[1]),10)[0::10])

plt.plot(data_1e_2[0][0::10], moving_average(data_1e_2[1]*0.2/(data_1e_2[2]-data_1e_2[1]),10)[0::10])
plt.plot(data_3e_2[0][0::10], moving_average(data_3e_2[1]*0.2/(data_3e_2[2]-data_3e_2[1]),10)[0::10] )
plt.legend( ['no ISAT',"other",r"$k = 1\times 10^{-3}$",r"$k = 3\times 10^{-3}$",r"$k = 1\times 10^{-2}$",r"$k = 3\times 10^{-2}$"],loc='upper left',frameon=False,fontsize=13)

plt.xlim(0,2e-6)
plt.ylim(0,3.5)
#plt.legend( ['VLE time', 'Others'],loc='upper left',frameon=False,fontsize=13)
plt.xlabel("physical time (s)", fontsize=12)
plt.ylabel("cpu time (s)", fontsize=12)

plt.savefig(folder + "//time_noISAT.png", dpi=400 , bbox_inches='tight')

