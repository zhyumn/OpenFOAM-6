#import Foam
#import VLE
from VTFlash import *
from utilities import *
from SpeciesData import *
import numpy as np
#import os

#print("inside the code in slurm")
# #set conditions v, U, Z, T, P
Species_names = ['C6H14','C12H26']
Species = [1,3] #species indices

#datarange
Yn2 = [0.6] #np.linspace(0.9,0.5,10)
T_space = np.linspace(560,600,50)
#U_space = np.linspace(-5650.61,-500,10)
vol_space = np.linspace(0.0008,0.0004,40)
Zm = [0,0]      #mole fractions of mixtures

#file name
f = open("pref_evap_c6_c12.txt", "w")
#f.write("T\tv\tZ1\tZ2\tP\tx1\tx2\ty1\ty2\n")

# #TP side
# ofpath = os.getenv('FOAM_DIR')
# mix = VLE.solver_new(ofpath+"/user/etc/caseDicts/VLEcase/")
# mix.reset_specie(["N2","C6H14"])

#mix.setX([0.5, 0.5])
#mix.setT(350)
#mix.setP(50e5)
#mix.setKinit(mix.K())
#mix.setTPn_flag(VLE.TPN_TPD_Tud)

# for Tcount in range(len(T_space)):
#     for volcount in range(len(vol_space)):
#         for ncount in range(len(Yn2)): 

for ncount in range(len(Yn2)): 
    for volcount in range(len(vol_space)):
        for Tcount in range(len(T_space)): 
            Zm[0] = Yn2[ncount]
            Zm[1] = 1.0 - Yn2[ncount] 
            Mavg = 0
            for i in range(len(Species)):
                Mavg += Zm[i]/SpeciesData[Species_names[i]]['W']
            Mavg = 1/Mavg
            Z = np.zeros(len(Species))
            for i in range(len(Species)):
            	Z[i] = 0.5 #Zm[i]*Mavg/SpeciesData[Species_names[i]]['W']
            
            print(Z)

            T = T_space[Tcount] 
            v = vol_space[volcount]   

            P = 20e5 
            R = 8.314

            x,y,vaporfrac,vl,vv = VTFlash(v,T,Z,P,Species,Species_names)

            P = P_PR(T,vv,y,Species,Species_names)

            print('T (K) ',T)
            print('v (m3/mole) ',v)
            print('Z (mole fraction)',Z)
            print('vaporfrac = '+str(vaporfrac))
            print('x = '+str(x))
            print('y = '+str(y))
            print('pressure - l (bar) = '+str(P_PR(T,vl,x,Species,Species_names)/100000))
            print('pressure - v (bar) = '+str(P_PR(T,vv,y,Species,Species_names)/100000))

            Pl = P_PR(T,vl,y,Species,Species_names)
            Pv = P_PR(T,vv,y,Species,Species_names)

            if(vaporfrac == 1):
                P = Pv #P_PR(T,vv,y,Species,Species_names)
            elif(vaporfrac == 0):
                P = Pl #P_PR(T,vl,y,Species,Species_names)
            else:
                P = Pv #P_PR(T,vv,y,Species,Species_names)
            
            if((np.isnan(Pl) and vaporfrac == 0) or (np.isnan(Pv) and vaporfrac == 1) or (P<=0)):
                print('nan')
            elif(vaporfrac > 0 and vaporfrac < 1 and (x[0]<=x[1])):
                f.write(str(T)+"\t"+str(v)+"\t"+str(Z[0])+"\t"+str(Z[1])+"\t"+str(P/100000)+"\t"+str(x[0])+"\t"+str(x[1])+"\t"+str(y[0])+"\t"+str(y[1])+"\t"+str(vaporfrac)+"\n")

            print("\n")

f.close()
