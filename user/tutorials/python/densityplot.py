
import Foam
import VLE
import matplotlib.pyplot as plt
import numpy as np
import os

ofpath = os.getenv('FOAM_DIR')
mix = VLE.solver_new(ofpath+"/user/etc/caseDicts/VLEcase/")
mix.reset_specie(["CO2","H2O"])

mix.setX([0.5, 0.5])
mix.setT(350)
mix.setP(50e5)
mix.setKinit(mix.K())
mix.setTPn_flag(VLE.TPN_TPD_Tud)

T = []
pressure = [5e5, 10e5,50e5,100e5,200e5]
rho = [ [] for i in range(len(pressure))]
cp = [ [] for i in range(len(pressure))]
vf = [ [] for i in range(len(pressure))]


for it in np.arange(100,600,1):
    mix.setT(float(it))
    T.append(it)
    for i in range(len(pressure)):
        mix.setP(pressure[i])
        vf[i].append(mix.vaporfra())
        cp[i].append(mix.Cp())
        rho[i].append(mix.rho())
fig, ax1 = plt.subplots()
for i in range(len(pressure)):
    ax1.plot(T,rho[i])

plt.legend(["5 bar","10 bar","50 bar","100 bar","200 bar"],frameon=False,labelspacing=0,fontsize=24,bbox_to_anchor=(0.77, 0.75),loc='center',handlelength=1.2)
plt.ylabel("Density $(kg/m^3)$", fontsize=24)
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
plt.xlabel("T (K)", fontsize=24)
ax2 = ax1.twinx()

for i in range(len(pressure)):
    ax2.plot(T,vf[i],"--")
plt.ylim(0,2)

plt.savefig("density_vf.png", dpi=300,bbox_inches='tight')

