import numpy as np
import os
import matplotlib.pyplot as plt


def get_number(s):
  s = s.replace("processor", "")
  return int(s)

def read(folder, name):
  proc_file = os.listdir(folder+"/CPU_data" + name)
  proc_file.sort(key=get_number)
  a = []

  fn = os.listdir(folder+"/CPU_data" + name+"/processor0")[0]

  for i in proc_file:
      file_name = folder+"/CPU_data" + name +"/"+ i +"/"+fn+"/"+fn+"/cpu_total.out"
      a.append([[float(s1) for s1 in s.split()] for s in open(file_name, 'r').read().split("\n")[:-1]])
  data = np.array(a)
  return data

folder = os.path.split(os.path.realpath(__file__))[0]

name = "64_np_8-8_1"

data_6_np = read(folder,"/"+name)

#data_4_1_nlb = read(folder,"/4_1_nlb" )
#s_4_1_nlb = [np.sum(i[:,1]) for i in data_4_1_nlb]
#s_4_1_nlb.append(np.sum(np.array(s_4_1_nlb)))

#data_4_1_lb = read(folder, "/4_1_lb_50")
#s_4_1_lb = [np.sum(i[:,1]) for i in data_4_1_lb]
#s_4_1_lb .append(np.sum(np.array(s_4_1_lb)))
'''
data_4_1_lb_new = read(folder, "/6_1_lb_50")
s_4_1_lb_new = [np.sum(i[:,1]) for i in data_4_1_lb_new]
s_4_1_lb_new .append(np.sum(np.array(s_4_1_lb_new)))

data_4_1_lb_new2= read(folder, "/6_1_lb_100")
s_4_1_lb_new2 = [np.sum(i[:,1]) for i in data_4_1_lb_new2]
s_4_1_lb_new2 .append(np.sum(np.array(s_4_1_lb_new2)))

data_4_1_lb_new3= read(folder, "/6_1_lb_100_2")
s_4_1_lb_new3 = [np.sum(i[:,1]) for i in data_4_1_lb_new3]
s_4_1_lb_new3 .append(np.sum(np.array(s_4_1_lb_new3)))
'''
'''
plt.plot(data_4_1_nlb[0,:,0],data_4_1_nlb[0,:,1])
plt.plot(data_4_1_nlb[1,:,0],data_4_1_nlb[1,:,1])
plt.plot(data_4_1_nlb[2,:,0],data_4_1_nlb[2,:,1])
plt.plot(data_4_1_nlb[3,:,0],data_4_1_nlb[3,:,1])
plt.plot(data_4_1_nlb[3,:,0],data_4_1_nlb[0,:,1]+data_4_1_nlb[1,:,1]+data_4_1_nlb[2,:,1]+data_4_1_nlb[3,:,1],"--")

plt.savefig(folder + "//time1.png", dpi=400 , bbox_inches='tight')
plt.cla()

plt.plot(data_4_1_lb[0,:,0],data_4_1_lb[0,:,1])
plt.plot(data_4_1_lb[1,:,0],data_4_1_lb[1,:,1])
plt.plot(data_4_1_lb[2,:,0],data_4_1_lb[2,:,1])
plt.plot(data_4_1_lb[3,:,0],data_4_1_lb[3,:,1])
plt.plot(data_4_1_nlb[3,:,0],data_4_1_nlb[3,:,1])
plt.plot(data_4_1_lb[3,:,0],data_4_1_lb[0,:,1]+data_4_1_lb[1,:,1]+data_4_1_lb[2,:,1]+data_4_1_lb[3,:,1],"--")
plt.plot(data_4_1_nlb[3,:,0],data_4_1_nlb[0,:,1]+data_4_1_nlb[1,:,1]+data_4_1_nlb[2,:,1]+data_4_1_nlb[3,:,1],"--")
plt.savefig(folder + "//time2.png", dpi=400 , bbox_inches='tight')

plt.cla()
'''



#index = np.arange(4)
#bar_width = 0.35
#plt.bar(index,s_4_1_nlb,bar_width,label="No load balancing")
#plt.bar(index+bar_width,s_4_1_lb,bar_width,label="load balancing" )
#plt.bar(index+2*bar_width,s_4_1_lb_new,bar_width )
#plt.bar(index+3*bar_width,s_4_1_lb_new2,bar_width )
#plt.bar(index+4*bar_width,s_4_1_lb_new3,bar_width )
#print(s_4_1_lb[4]/s_4_1_nlb[4])
#plt.legend()
#plt.xticks(index, ["0","1","2","3"])
#plt.xlabel("MPI rank ", fontsize=12)
#plt.ylabel("cpu time (s)", fontsize=12)
#plt.savefig(folder + "//time5.png", dpi=400 , bbox_inches='tight')

data_6_np_max = np.array([np.max(data_6_np[:,i,1]) for i in range(len(data_6_np[0]))])


print(np.sum(data_6_np_max),np.sum(data_6_np[10,:,1]))