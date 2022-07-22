
import Foam
import VLE
import csv
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from VLE import DoubleVector as vec
import os

import sys

#from user.run.Python.TP_diagram import P
print(sys.path)
ofpath = os.getenv('FOAM_DIR')

mix = VLE.solver_new(ofpath+"/user/etc/caseDicts/VLEcase/")

mix.reset_specie(["CO2","H2O","CH4" ,"O2"])
mix.setX([0.9, 0.01,0.045,0.045])

mix.setTPn_flag(VLE.TPN_TPD)
mix.setTPn_flag(VLE.TPN_old)
T = []
rho0 = []
rho1 = []
rho2 = []
rho3 = []
rho4 = []
cp0 =[]
cp1 =[]
cp2 =[]
cp3 =[]
cp4 =[]
vf0=[]
vf1=[]
vf2=[]
vf3=[]
vf4=[]

mix.setTPn_flag(VLE.TPN_TPD_Tud)

for it in np.arange(100,600,1):
    #print(it)
    mix.setT(float(it))
    T.append(it)
    mix.setP(5e5)
    vf0.append(mix.vaporfra())
    cp0.append(mix.Cp())
    rho0.append(mix.rho())
    #print(5e5)
    
    mix.setP(10e5)

    vf1.append(mix.vaporfra())
    cp1.append(mix.Cp())
    rho1.append(mix.rho())
    #print(10e5)
    mix.setP(50e5)

    vf2.append(mix.vaporfra())
    cp2.append(mix.Cp())
    rho2.append(mix.rho())
    #print(50e5)
    mix.setP(100e5)

    vf3.append(mix.vaporfra())
    cp3.append(mix.Cp())
    rho3.append(mix.rho())
    #print(100e5)
    mix.setP(200e5)
    vf4.append(mix.vaporfra())
    rho4.append(mix.rho())
    #print(200e5)
fig, ax1 = plt.subplots()
ax1.plot(T,rho0)
ax1.plot(T,rho1)
ax1.plot(T,rho2)
ax1.plot(T,rho3)
ax1.plot(T,rho4)
plt.legend(["5 bar","10 bar","50 bar","100 bar","200 bar"],frameon=False,labelspacing=0,fontsize=24,bbox_to_anchor=(0.77, 0.75),loc='center',handlelength=1.2)
plt.ylabel("Density $(kg/m^3)$", fontsize=24)
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
plt.xlabel("T (K)", fontsize=24)
ax2 = ax1.twinx()
ax2.plot(T,vf0,"--")
ax2.plot(T,vf1,"--")
ax2.plot(T,vf2,'--')

plt.ylim(0,2)
#plt.legend(["old","new"],frameon=False,labelspacing=0,fontsize=24,bbox_to_anchor=(0.77, 0.75),loc='center',handlelength=1.2)
pair = [(T[i],vf1[i]) for i in range(len(T))]
#print(pair)
plt.savefig("vfdbug.png", dpi=300,bbox_inches='tight')
"""

mix.reset_specie(["CO2","H2O","CH4" ,"O2"])
mix.setX([0.9, 0.01,0.045,0.045])
mix.setP(5e5)
mix.setT(109)

mix.setTPn_flag(VLE.TPN_old)
print(mix.vaporfra())
mix.setTPn_flag(VLE.TPN_TPD_Tud)
print(mix.vaporfra())
#mix.Cp()
"""
