import numpy as np
from math import pi
from os.path import exists
import os

with open("Result20220727T222437722Z.txt", "r") as file:
    rawdata=list(file)
    data=np.array([[float(s) for s in l.split(';')] for l in rawdata[1:]])

data[:,0]/=1000
r=4.5e-5
rho=643.25

data[:,1]/=1000*rho*pi*r*r

file_exists = exists("0/U.table")
if file_exists:
    os.system("mv 0/U.table 0/U.table.back")

stringout = ["("+"%.7f"% l[0] + " (" +"%.7f"% l[1]+" 0 0))\n"  for l in data]
with open("0/U.table", "w") as file:
    for i in range(7501):
        file.write(stringout[i])



