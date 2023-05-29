import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("pref_evap_c6_c12.txt",delimiter='\t',dtype=float)

T = data[:,0]
P = data[:,4]

x1 = data[:,5]
x2 = data[:,6]
y1 = data[:,7]
y2 = data[:,8]
vf = data[:,9]

plt.scatter(T,P,c=vf,marker='s',cmap='plasma')
#plt.scatter(T,y2,c=vf,marker='^',cmap='plasma')
#plt.title("Preferential Evaporation (C6H14 - C12H24)",fontsize=15)
# plt.legend(['Liquid Phase','Gas Phase'],fontsize=12)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('Temperature (K)',font="Helvetica",fontsize=15)
plt.ylabel('Pressure (bar)',font="Helvetica",fontsize=15)
plt.colorbar()
plt.savefig("pref-evap_P.png", dpi=300,bbox_inches='tight')