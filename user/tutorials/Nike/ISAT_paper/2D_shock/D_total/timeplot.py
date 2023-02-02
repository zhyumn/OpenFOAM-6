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
data_noISAT=read_data(folder + "/time_noISAT")
data_1e_2 =read_data(folder + "/time_1e_2")

file = open(folder+ "/time_noISAT/VLEtime", 'r')
data = [s.split(",") for s in file.read().split("\n")[:-1]]
data = np.array([[float(s2) for s2 in s1] for s1 in data])

file = open(folder+ "/time_noISAT/cputotal", 'r')
data2 = [s.split(",") for s in file.read().split("\n")[:-1]]
data2 = np.array([[float(s2) for s2 in s1] for s1 in data2])

folder = os.path.split(os.path.realpath(__file__))[0]
file = open(folder+ "/time_1e_3/VLEtime", 'r')
data3 = [s.split(",") for s in file.read().split("\n")[:-1]]
data3 = np.array([[float(s2) for s2 in s1] for s1 in data3])

folder = os.path.split(os.path.realpath(__file__))[0]
file = open(folder+ "/time_1e_2/VLEtime", 'r')
data4 = np.array([[float(s2) for s2 in s1] for s1 in [s.split(",") for s in file.read().split("\n")[:-1]]])

folder = os.path.split(os.path.realpath(__file__))[0]
file = open(folder+ "/time_1e_4/VLEtime", 'r')
data5 = np.array([[float(s2) for s2 in s1] for s1 in [s.split(",") for s in file.read().split("\n")[:-1]]])

folder = os.path.split(os.path.realpath(__file__))[0]
file = open(folder+ "/time_3e_2/VLEtime", 'r')
data6 = np.array([[float(s2) for s2 in s1] for s1 in [s.split(",") for s in file.read().split("\n")[:-1]]])

folder = os.path.split(os.path.realpath(__file__))[0]
file = open(folder+ "/time_3e_3/VLEtime", 'r')
data7 = np.array([[float(s2) for s2 in s1] for s1 in [s.split(",") for s in file.read().split("\n")[:-1]]])



def moving_average(interval,windowsize):
    window = np.ones(int(windowsize))/float(windowsize)
    re = np.convolve(interval,window,'same')
    return re
#----------------------------------------------
plt.plot(data_noISAT[0][0::10], moving_average(data_noISAT[1],10)[0::10],  color= "r")
plt.plot(data3[:,0][0::10], moving_average(data3[:,1],10)[0::10]  )
plt.plot(data7[:,0][0::10], moving_average(data7[:,1],10)[0::10]  ,  color= "black")

#plt.plot(data4[:,0][0::10], moving_average(data4[:,1],10)[0::10]  )
plt.plot(data_1e_2[0][0::10], moving_average(data_1e_2[1],10)[0::10]  )
plt.plot(data6[:,0][0::10], moving_average(data6[:,1],10)[0::10]  )


#plt.plot(data5[:,0][0::10], moving_average(data5[:,1],30)[0::10]  )



plt.plot(data_noISAT[0][0::10], moving_average(data_noISAT[2]-data_noISAT[1],10)[0::10],'--', color= "b")


plt.xlim(0,1.5e-6)
plt.ylim(0,3.5)
plt.legend( ['no ISAT',r"$k = 1\times 10^{-3}$",r"$k = 3\times 10^{-3}$",r"$k = 1\times 10^{-2}$",r"$k = 3\times 10^{-2}$"],loc='upper left',frameon=False,fontsize=13)
plt.xlabel("physical time (s)", fontsize=12)
plt.ylabel("cpu time (s)", fontsize=12)
print("fig saved")
plt.savefig(folder + "//time.png", dpi=400 , bbox_inches='tight')
#----------------------------------------------
plt.cla()

plt.plot(data[:,0][0::10], moving_average(data[:,1],10)[0::10],  color= "r")
plt.plot(data2[:,0][0::10], moving_average(data2[:,1]-data[:,1],10)[0::10],'--', color= "b")

plt.xlim(0,1.5e-6)
plt.ylim(0,3.5)
plt.legend( ['VLE time', 'Others'],loc='upper left',frameon=False,fontsize=13)
plt.xlabel("physical time (s)", fontsize=12)
plt.ylabel("cpu time (s)", fontsize=12)

plt.savefig(folder + "//time_noISAT.png", dpi=400 , bbox_inches='tight')

#----------------------------------------------
NN=1954

OT=np.sum(data2[:NN,1])-np.sum(data[:NN,1])
T1=np.sum(data[:NN,1])
T2=np.sum(data4[:NN,1])
T3=np.sum(data3[:NN,1])
T4=np.sum(data5[:NN,1])
T5=np.sum(data6[:NN,1])
T6=np.sum(data7[:NN,1])
#----------------------------------------------

plt.cla()
plt.plot([1e-3,3e-3,1e-2,3e-2], [T5,T2,T6,T3])
plt.plot([1e-3,3e-3,1e-2,3e-2], [OT+T5,OT+T2,OT+T6,OT+T3])
#plt.plot([1e-3,3e-3,1e-2,3e-2], [T5,T2,T6,T3])
plt.savefig(folder + "//t.png", dpi=400 , bbox_inches='tight')
#----------------------------------------------
plt.cla()

plt.plot(data[:,0],moving_average(data[:2646,1]/data6[:2646,1],10))

plt.xlim(1.4e-6,1.5e-6)
plt.ylim(0,8)
plt.savefig(folder + "//ratio.png", dpi=400 , bbox_inches='tight')
#----------------------------------------------
plt.cla()

plt.pie([np.sum(data[:,1]),np.sum(data2[:,1]-data[:,1])],labels=["VLE","Others"],autopct='%1.1f%%',startangle=90,textprops={'fontsize': 14})
plt.savefig(folder + "//timepercent.png", dpi=400 , bbox_inches='tight')


#----------------------------------------------


print(data[:,1]/data7[:,1])
print("-----------")
print(T1)
print(T5)
print(T2)
print(T6)
print(T3)
print(T4)

print("-----------")
print(T1/T5)
print(T1/T2)
print(T1/T6)
print(T1/T3)
print(T1/T4)

print("-----------")

print(OT+T1)
print(OT+T5)
print(OT+T2)
print(OT+T6)
print(OT+T3)
print(OT+T4)

print("-----------")
print((OT+T1)/(OT+T5))
print((OT+T1)/(OT+T2))
print((OT+T1)/(OT+T6))
print((OT+T1)/(OT+T3))
print((OT+T1)/(OT+T4))

print("-----------")

"""


tT=1e0 
line_T_max, = plt.plot(data[:,0], data[:,1]/tT, "--",color= "b")
line_T_mean, = plt.plot(data2[:,0], data2[:,1]/tT,color= "b")
plt.yscale('log')


te=1e2

folder = os.path.split(os.path.realpath(__file__))[0]
file = open(folder+ "/postProcessing/cellMax(mag(subtract(e_ex,e)))/0/volFieldValue.dat", 'r')
data = [s.split("\t") for s in file.read().split("\n")[6:-1]]
data = np.array([[float(s2) for s2 in s1] for s1 in data])

file = open(folder+ "/postProcessing/volAverage(mag(subtract(e_ex,e)))/0/volFieldValue.dat", 'r')
data2 = [s.split("\t") for s in file.read().split("\n")[6:-1]]
data2 = np.array([[float(s2) for s2 in s1] for s1 in data2])

line_e_max, =plt.plot(data[:,0], data[:,1]/te, "--", color= "r")
line_e_mean, =plt.plot(data2[:,0], data2[:,1]/te,color= "r")
plt.legend([line_T_mean, line_e_mean], ['Temperature', 'Internal energy'],loc='lower right',frameon=False,fontsize=13)
plt.xlabel("time/s", fontsize=12)
plt.axhline(y = 1, color='black', linestyle='dotted')
plt.text(0.1e-6,5,"Maximum error", fontsize=14)
plt.text(0.5e-6,2e-4,"Average error", fontsize=14)
print("fig saved")
plt.savefig(folder + "//error.png", dpi=400 , bbox_inches='tight')
"""