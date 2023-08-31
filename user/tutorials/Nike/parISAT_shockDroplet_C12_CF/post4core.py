import numpy as np
import os
import matplotlib.pyplot as plt
folder = os.path.split(os.path.realpath(__file__))[0]


def read_data(f, n):

    ret = []
    file = open(f+"/processor" + str(0)+"/ISAT_VLE/ISAT_VLE/VLEtime", 'r')
    data_t = [s.split(",") for s in file.read().split("\n")[:-1]]
    data_t = np.array([[float(s2) for s2 in s1] for s1 in data_t])
    ret.append(data_t[:, 0])

    for i in range(n):
        file = open(f+"/processor" + str(i)+"/ISAT_VLE/ISAT_VLE/VLEtime", 'r')
        data_t = [s.split(",") for s in file.read().split("\n")[:-1]]
        data_t = np.array([[float(s2) for s2 in s1] for s1 in data_t])
        ret.append(data_t[:, 1])

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
data4 = read_data(folder+"/CPU_data/np_4_4-1", 4)
data5 = read_data(folder+"/CPU_data/nlb_4_4-1", 4)
data6 = read_data(folder+"/CPU_data/lb_4_4-1_new", 4)
data7 = read_data(folder+"/CPU_data/lb_4_4-1_new9", 4)

data8 = read_data(folder+"/CPU_data/lb_4_2-2", 4)
data9 = read_data(folder+"/CPU_data/nlb_4_2-2", 4)

data10 = read_data(folder+"/CPU_data/np_4_2-2", 4)

# ----------------------------------------------


def moving_average(interval, windowsize):
    window = np.ones(int(windowsize))/float(windowsize)
    re = np.convolve(interval, window, 'same')
    return re
def clean_data(x,y):
    for i in range(len(y)):
        if y[i]<=0:
            y[i]=y[i-1]
            x[i]=x[i-1]
    return x,y
    


# ----------------------------------------------
'''
file = open(folder +"/ptime_data", 'w')
file.write("\n".join([",".join([str(s1) for s1 in s2]) for s2 in data.transpose()]))


plt.plot(data[0][0::10], moving_average(data[1],40)[0::10],  color= "b")
#plt.plot(data1[0][0::10], moving_average(data1[2]-data1[1],10)[0::10],  color= "r")

plt.savefig(folder + "//time_noISAT.png", dpi=400 , bbox_inches='tight') '''
plt.xlim(0, 1e-6)
''' plt.plot(data[0][0::10], moving_average(data[1], 40)[0::10],  color="r")
plt.plot(data[0][0::10], moving_average(data[2], 40)[0::10],  color="r")

plt.plot(data2[0][0::10], moving_average(data2[1], 40)[0::10],  color="b")
plt.plot(data2[0][0::10], moving_average(data2[2], 40)[0::10],  color="b") '''

plt.plot(data3[0][0::10], moving_average(data3[1], 40)[0::10],  color="black")
plt.plot(data3[0][0::10], moving_average(data3[2], 40)[0::10],  color="black")

plt.plot(data4[0][0::10], moving_average(data4[1], 40)[0::10], "--", color="b")
plt.plot(data4[0][0::10], moving_average(data4[2], 40)[0::10], "--", color="b")
plt.plot(data4[0][0::10], moving_average(data4[3], 40)[0::10], "--", color="b")
plt.plot(data4[0][0::10], moving_average(data4[4], 40)[0::10], "--", color="b")

plt.plot(data5[0][0::10], moving_average(data5[1], 40)[0::10], "--", color="r")
plt.plot(data5[0][0::10], moving_average(data5[2], 40)[0::10], "--", color="r")
plt.plot(data5[0][0::10], moving_average(data5[3], 40)[0::10], "--", color="r")
plt.plot(data5[0][0::10], moving_average(data5[4], 40)[0::10], "--", color="r") 

''' plt.plot(data6[0][0::10], moving_average(data6[1], 40)[0::10], "--", color="black")
plt.plot(data6[0][0::10], moving_average(data6[2], 40)[0::10], "--", color="black")
plt.plot(data6[0][0::10], moving_average(data6[3], 40)[0::10], "--", color="black")
plt.plot(data6[0][0::10], moving_average(data6[4], 40)[0::10], "--", color="black") '''

x,y=clean_data(data7[0],data7[4])
#plt.plot(data7[0][0::10], moving_average(data7[4], 40)[0::10], "--", color="green")
plt.plot(x[0::10], moving_average(y, 40)[0::10], "--", color="black")

plt.plot(data8[0][0::10], moving_average(data8[3], 40)[0::10], "--", color="green")

plt.plot(data9[0][0::10], moving_average(data9[1], 40)[0::10], "--", color="y")
plt.plot(data9[0][0::10], moving_average(data9[2], 40)[0::10], "--", color="y")
plt.plot(data9[0][0::10], moving_average(data9[3], 40)[0::10], "--", color="y")
plt.plot(data9[0][0::10], moving_average(data9[4], 40)[0::10], "--", color="y") 

plt.plot(data10[0][0::10], moving_average(data10[1], 40)[0::10], "--", color="green")
plt.plot(data10[0][0::10], moving_average(data10[2], 40)[0::10], "--", color="green")
plt.plot(data10[0][0::10], moving_average(data10[3], 40)[0::10], "--", color="green")
plt.plot(data10[0][0::10], moving_average(data10[4], 40)[0::10], "--", color="green") 
plt.savefig(folder + "//time_noISAT.png", dpi=400, bbox_inches='tight')
