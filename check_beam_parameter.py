import math
import numpy as np
import matplotlib.pyplot as plt


Theta = 0.5*120/1500 #rad
lambdaL = 0.8 #um
w0 =  10 #um

#solve for M^2 in experimental parameter
M2 = Theta * math.pi * w0/ lambdaL
print("M2:", M2)



ry = math.pi * (w0 ** 2) / (lambdaL)
print(ry, "ry in um without M2")
ry = math.pi * (w0**2)*M2/lambdaL
print(ry, "ry in um with *M2_theo")
M2_ex = 2.2
ry = math.pi * (w0**2)*M2_ex/lambdaL
print(ry, "ry in um with *M2_exp")

diameter = 120
focalLength = 1500

waist0 = (2 * lambdaL / math.pi) * focalLength / diameter
print("w0 theo without M^2 by divergence", waist0)

a0= ((6*1E19*(lambdaL)**2 )/(1.38*1E18))**0.5
print("a0 for 6x10^19W/cm12", a0)