import numpy as np
import os
import matplotlib.pyplot as plt
folder = os.path.split(os.path.realpath(__file__))[0]


def read_data(f, n):

    ret = []
    file = open(f+"/processor" + str(0)+"/ISAT_VLE/ISAT_VLE/VLEtime", 'r')
    data_t = [s.split(",") for s in file.read().split("\n")[:-1]]
    data_t = np.array([[float(s2) for s2 in s1] for s1 in data_t])
    ret.append(data_t[:,0])
    
    for i in range(n):
        file = open(f+"/processor" + str(i)+"/ISAT_VLE/ISAT_VLE/VLEtime", 'r')
        data_t = [s.split(",") for s in file.read().split("\n")[:-1]]
        data_t = np.array([[float(s2) for s2 in s1] for s1 in data_t])
        ret.append(data_t[:,1])

    file = open(f+"/processor" + str(0) + "/ISAT_VLE/ISAT_VLE/cputotal", 'r')
    data2 = [s.split(",") for s in file.read().split("\n")[:-1]]
    data2 = np.array([[float(s2) for s2 in s1] for s1 in data2])
    ret.append(data2[:, 1])
    ret = np.array(ret)
    return ret
    return np.array([data[:, 0], data[:, 1], data2[:, 1]])


data = read_data(folder+"/CPU_data/nlb_2_new", 2)
data2 = read_data(folder+"/CPU_data/np_2", 2)
data3 = read_data(folder+"/CPU_data/lb_2", 2)


 #----------------------------------------------
def moving_average(interval,windowsize):
    window = np.ones(int(windowsize))/float(windowsize)
    re = np.convolve(interval,window,'same')
    return re
#----------------------------------------------
'''
file = open(folder +"/ptime_data", 'w')
file.write("\n".join([",".join([str(s1) for s1 in s2]) for s2 in data.transpose()]))


plt.plot(data[0][0::10], moving_average(data[1],40)[0::10],  color= "b")
#plt.plot(data1[0][0::10], moving_average(data1[2]-data1[1],10)[0::10],  color= "r")

plt.savefig(folder + "//time_noISAT.png", dpi=400 , bbox_inches='tight') '''
plt.xlim(0,1e-6)
plt.plot(data[0][0::10], moving_average(data[1],40)[0::10],  color= "r")
plt.plot(data[0][0::10], moving_average(data[2],40)[0::10],  color= "r")

plt.plot(data2[0][0::10], moving_average(data2[1],40)[0::10],  color= "b")
plt.plot(data2[0][0::10], moving_average(data2[2],40)[0::10],  color= "b")

plt.plot(data3[0][0::10], moving_average(data3[1],40)[0::10],  color= "black")
plt.plot(data3[0][0::10], moving_average(data3[2],40)[0::10],  color= "black")
plt.savefig(folder + "//time_noISAT.png", dpi=400 , bbox_inches='tight')