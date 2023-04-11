import numpy as np
import os
import matplotlib.pyplot as plt


# 定义一个提取数字的函数
def get_number(s):
  # 去掉字符串中的"process"部分
  s = s.replace("processor", "")
  # 将剩下的部分转换为整数
  return int(s)

def read(folder, name):
  proc_file = os.listdir(folder+"/CPU_data" + name)
  proc_file.sort(key=get_number)
  a = []

  for i in proc_file:
      file_name = folder+"/CPU_data" + name +"/"+ i +"/STDAC/STDAC/cpu_total.out"
      a.append([[float(s1) for s1 in s.split()] for s in open(file_name, 'r').read().split("\n")[:-1]])
  data = np.array(a)
  return data

folder = os.path.split(os.path.realpath(__file__))[0]



data_4_1_nlb = read(folder,"/6_1_nlb" )
s_4_1_nlb = [np.sum(i[:,1]) for i in data_4_1_nlb]
s_4_1_nlb.append(np.sum(np.array(s_4_1_nlb)))

data_4_1_lb = read(folder, "/6_1_lb_50_2")
s_4_1_lb = [np.sum(i[:,1]) for i in data_4_1_lb]
s_4_1_lb .append(np.sum(np.array(s_4_1_lb)))

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
index = np.arange(7)
bar_width = 0.1
plt.bar(index,s_4_1_nlb,bar_width)
plt.bar(index+bar_width,s_4_1_lb,bar_width )
plt.bar(index+2*bar_width,s_4_1_lb_new,bar_width )
plt.bar(index+3*bar_width,s_4_1_lb_new2,bar_width )
plt.bar(index+4*bar_width,s_4_1_lb_new3,bar_width )
print(s_4_1_lb[4]/s_4_1_nlb[4])
plt.savefig(folder + "//time4.png", dpi=400 , bbox_inches='tight')