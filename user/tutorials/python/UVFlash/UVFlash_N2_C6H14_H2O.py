from VTFlash import *
from utilities import *
import numpy as np

#set conditions v, U, Z, T, P
Z = [0.5,0.35,0.15]      #mole fractions of mixtures
Species_names = ['N2','C12H16']
Species = [0,3] #species indices from Species data to ensure correct binary interaction coefficient is being used
T = 400
v = 0.00051433 #V_PR(T,P,Z,Species,Species_names)       #specific volume of mixture v/moles
P = 40e5 #P_PR(T,v,Z,Species,Species_names) #1000000       #pressure
print('pr vol ',V_PR(T,P,Z,Species,Species_names))
R = 8.314

# u = Umix(v,T,P,Z,Species)
u = -58029.79719905442 #set this
cv = 141.25241198730748

error = 100

# #step 0 - find mixture T assuming single phase
T = T_SinglePhaseGuess(u,Z,P,v,Species,Species_names)

count = 0

while error > 1E-8:
    #step 1 - run vT flash using above T
    x,y,vaporfrac,vl,vv = VTFlash(v,T,Z,P,Species,Species_names)
    #print('vf ',vaporfrac)

    #step 2 - get umix and cv mix
    u_mix, cv_mix = u_cv_mix_real(x,y,vaporfrac,T,P,vl,vv,Species,Species_names)
    
    #step 3 - update T
    T = T + 0.1*(u - u_mix)/cv_mix

    #update pressure 
    P = P_PR(T,vv,y,Species,Species_names)

    #step 4 - estimate error
    error = np.abs((u-u_mix)/u)
    count += 1
    
    #print(count,u_mix,cv_mix,T)
    #print(error)

    if(count == 10000): #or vaporfrac == 1 or vaporfrac == 0):
        break

print('T (K) ',T)
print('v (m3/mole) ',v)
print('Z (mole fraction)',Z)
print('vaporfrac = '+str(vaporfrac))
print('x = '+str(x))
print('y = '+str(y))
print('pressure - l = '+str(P_PR(T,vl,x,Species,Species_names)/100000))
print('pressure - v = '+str(P_PR(T,vv,y,Species,Species_names)/100000))

