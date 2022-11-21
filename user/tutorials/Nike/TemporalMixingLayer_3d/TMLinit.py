# %%
import math
import matplotlib.pyplot as plt
import Foam
import VLE
import csv
# from matplotlib import cm
# from mpl_toolkits.mplot3d import Axes3D

import numpy as np
from VLE import DoubleVector as vec

import sys

# from user.run.Python.TP_diagram import P
print(sys.path)

mix1 = VLE.solver_new(
    "/home/zhy/new_drive/OpenFOAM-6/user/etc/caseDicts/VLEcase/")
mix2 = VLE.solver_new(
    "/home/zhy/new_drive/OpenFOAM-6/user/etc/caseDicts/VLEcase/")

mix1.reset_specie(["C6H14", "N2"])
mix2.reset_specie(["C6H14", "N2"])
mix1.setP(5e6)
mix2.setP(5e6)
C=np.array([0.9999, 0.0001])
N=np.array([0.0001, 0.9999])
mix1.setY(C)
mix2.setY(N)
mix1.setT(595)
mix2.setT(293)
mix1.setTPn_flag(VLE.TPN_TPD_Tud)
mix2.setTPn_flag(VLE.TPN_TPD_Tud)


rho1=mix1.rho()
rho2=mix2.rho()
rhom=0.5*(rho1+rho2)

mu1=mix1.mu()
print(mu1)
mu2=mix2.mu()
print(mu2)
mum=0.5*(mu1+mu2)

c1 = mix1.c()
c2 = mix2.c()
print(c1,c2)

w1= mix1.W()
w2= mix2.W()
print(w1,w2)

R =8.31446261815324
Mc=0.4

coef= math.sqrt(mix1.P()/mix2.P()*w1/w2*mix2.T()/mix1.T())
U1 = 2*Mc*c1/(1+c1/c2*coef)
U2 =-coef*U1
print(U1,U2)



dU = U1-U2
L1= 1.735e-5


delta = L1/29.16

Re0 = rhom*dU*delta/mum
print(Re0)


