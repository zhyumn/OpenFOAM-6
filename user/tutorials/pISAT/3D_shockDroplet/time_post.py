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

data = read_data(folder+"/CPU_data/lb_8_p2", 8)
#data2 = read_data(folder+"/CPU_data/nlb_8", 8)
data2 = read_data(folder+"/CPU_data/lb3_8_p", 8)
data3 = read_data(folder+"/CPU_data/np_8_new", 8)
data4 = read_data(folder+"/CPU_data/o_8", 8)

 #----------------------------------------------
def moving_average(interval,windowsize):
    interval = np.concatenate((interval,np.ones(int(windowsize))*interval[-1]),axis=None)
    window = np.ones(int(windowsize))/float(windowsize)
    re = np.convolve(interval,window,'same')
    return re[:-windowsize]
#----------------------------------------------

data_max = np.maximum.reduce(data[1:9])
data2_max = np.maximum.reduce(data2[1:9])
data3_max = np.maximum.reduce(data3[1:9])
data4_max = np.maximum.reduce(data4[1:9])

plt.plot(data[0][0::10]*1e6, moving_average(data[1],40)[0::10],  color= "black")
plt.plot(data2[0][0::10]*1e6, moving_average(data2[1],40)[0::10],  color= "b")
plt.plot(data3[0][0::10]*1e6, moving_average(data3_max,40)[0::10],  color= "red")
plt.plot(data4[0][0::10]*1e6, moving_average(data4_max,40)[0::10],  color= "g")

""" plt.plot(data3[0][0::10]*1e6, moving_average(data3[1],40)[0::10],  color= "red")
plt.plot(data3[0][0::10]*1e6, moving_average(data3[2],40)[0::10],  color= "red")
plt.plot(data3[0][0::10]*1e6, moving_average(data3[3],40)[0::10],  color= "red")
plt.plot(data3[0][0::10]*1e6, moving_average(data3[4],40)[0::10],  color= "red")
plt.plot(data3[0][0::10]*1e6, moving_average(data3[5],40)[0::10],  color= "red")
plt.plot(data3[0][0::10]*1e6, moving_average(data3[6],40)[0::10],  color= "red")
plt.plot(data3[0][0::10]*1e6, moving_average(data3[7],40)[0::10],  color= "red")
 """
plt.savefig(folder + "//time_test.png", dpi=400 , bbox_inches='tight')




print(np.sum(data_max))
print(np.sum(data2_max))
print(np.sum(data3_max))
print(np.sum(data4_max))