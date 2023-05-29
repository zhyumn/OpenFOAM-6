import Foam
import VLE
from VTFlash import *
from utilities import *
from SpeciesData import *
import numpy as np
import os

# #set conditions v, U, Z, T, P
Species_names = ['N2','C6H14']
Species = [0,1] #species indices

#datarange
Yn2 = np.linspace(0.9,0.3,20)
U_space = np.linspace(5986.526315789473,-4600,10)
vol_space = np.linspace(0.00116024,0.000311663,20)
Zm = [0,0]      #mole fractions of mixtures

#file name
f = open("data_gen_IC_shock_new.txt", "a")
f.write("U\tv\tZ1\tZ2\tT\tP\tvf\tc\n")

#TP side
ofpath = os.getenv('FOAM_DIR')
mix = VLE.solver_new(ofpath+"/user/etc/caseDicts/VLEcase/")
mix.reset_specie(["N2","C6H14"])

#mix.setX([0.5, 0.5])
#mix.setT(350)
#mix.setP(50e5)
#mix.setKinit(mix.K())
#mix.setTPn_flag(VLE.TPN_TPD_Tud)


for ucount in range(len(U_space)):
    for volcount in range(len(vol_space)):
        for ncount in range(len(Yn2)):

            Zm[0] = Yn2[ncount]
            Zm[1] = 1.0 - Yn2[ncount] 
            Mavg = 0
            for i in range(len(Species)):
                Mavg += Zm[i]/SpeciesData[Species_names[i]]['W']
            Mavg = 1/Mavg
            Z = np.zeros(len(Species))
            for i in range(len(Species)):
            	Z[i] = Zm[i]*Mavg/SpeciesData[Species_names[i]]['W']
            
            print(Z)
            T = 400

            v = vol_space[volcount]   

            P = 27e5 
            R = 8.314

            
            u = U_space[ucount]
            cv = 141.25241198730748

            error = 100

            # #step 0 - find mixture T assuming single phase
            T = T_SinglePhaseGuess(u,Z,P,v,Species,Species_names)
            print('T guess:\t',T)

            count = 0
            #using numerical CV
            T_old_1 = T 
            x,y,vaporfrac,vl,vv = VTFlash(v,T_old_1,Z,P,Species,Species_names)
            u_mix_1, cv_mix_old = u_cv_mix_real(x,y,vaporfrac,T_old_1,P,vl,vv,Species,Species_names)

            #finding u when T+1
            T_old_2 = T+1
            x,y,vaporfrac,vl,vv = VTFlash(v,T_old_2,Z,P,Species,Species_names)
            u_mix_2, cv_mix_old = u_cv_mix_real(x,y,vaporfrac,T_old_2,P,vl,vv,Species,Species_names)

            u_mix_old_1 = u_mix_1
            u_mix_old_2 = u_mix_2

            cv_mix = (u_mix_old_2 - u_mix_old_1)/(T_old_2 - T_old_1)
            error_old = 100

            alpha = 0.01

            while error > 1E-4:
                #step 1 - run vT flash using above T
                x,y,vaporfrac,vl,vv = VTFlash(v,T,Z,P,Species,Species_names)
                #print('vf ',vaporfrac)

                #step 2 - get umix and cv mix
                u_mix, cv_mix = u_cv_mix_real(x,y,vaporfrac,T,P,vl,vv,Species,Species_names)
                dT = 1
                u_mix_dT, cv_mix_dT = u_cv_mix_real(x,y,vaporfrac,T+dT,P,vl,vv,Species,Species_names)
                cv_mix_num = np.abs((u_mix_dT - u_mix)/dT) #np.abs((u_mix_old_2 - u_mix_old_1)/(T_old_2 - T_old_1 + 1E-20))
                
                #step 3 - update T
                T = T + alpha*(u - u_mix)/cv_mix_num

                #update pressure 
                P = P_PR(T,vv,y,Species,Species_names)

                #step 4 - estimate error
                error = np.abs((u-u_mix)/u)
                count += 1

                #old data saved
                T_old_1 = T_old_2 
                T_old_2 = T
                u_mix_old_1 = u_mix_old_2
                u_mix_old_2 = u_mix
                #print(count,u_mix,cv_mix,cv_mix_num,T,error)
                #print(error)

                # if(np.abs(T_old_1 - T_old_2)<1E-5 or np.abs(error - error_old)<1E-4):
                #     break
                
                error_old = error

                if(count == 1000):
                    alpha = alpha/10

                if(count == 2000): #or vaporfrac == 1 or vaporfrac == 0):
                    break


            print('u ',u)
            print('v (m3/mole) ',v)
            print('Z (mole fraction)',Z)
            print('T (K) ',T)
            print('vaporfrac = '+str(vaporfrac))
            print('x = '+str(x))
            print('y = '+str(y))
            print('pressure - l (bar) = '+str(P_PR(T,vl,x,Species,Species_names)/100000))
            print('pressure - v (bar) = '+str(P_PR(T,vv,y,Species,Species_names)/100000))

            if(vaporfrac == 1):
                P = P_PR(T,vv,y,Species,Species_names)
            elif(vaporfrac == 0):
                P = P_PR(T,vl,y,Species,Species_names)
            else:
                P = P_PR(T,vv,y,Species,Species_names)
           
            mix.setX([Z[0],Z[1]])
            mix.setT(T)
            mix.setP(P)
            mix.setKinit(mix.K())
            mix.setTPn_flag(VLE.TPN_TPD_Tud)		
            c = mix.c()

            print('speed of sound = '+str(c))

            f.write(str(u)+"\t"+str(v)+"\t"+str(Z[0])+"\t"+str(Z[1])+"\t"+str(T)+"\t"+str(P)+"\t"+str(vaporfrac)+"\t"+str(c)+"\n")

            print("\n")
            f.flush()

f.close()
