#============================ Module Importation ==================================================
import collections as co
import random
import statistics

from scipy import optimize
import matplotlib.pyplot as plt
import numpy as np

import pickle
from copy import deepcopy

import os
import sys
import subprocess
import time
#==================================================================================================
#============================ Global Definition ===================================================
# TD=[1e-6, 2e-6, 3e-6, 4e-6, 5e-6, 6e-6, 7e-6, 8e-6, 9e-6] # ----Throat Parameter is Diameter here
# TS=[0, 3, 4, 5, 6] # -------------------------------------------Throat Crosssection Polygon Shape
Throat=co.namedtuple('T', 'D, S') # ------------------Throat Type: Diameter, Shape of Cross Section
ThroatChoice={0:{1:Throat(1e-6,0), 2:Throat(2e-6,0), 3:Throat(3e-6,0), 4:Throat(4e-6,0), 5:Throat(5e-6,0), 6:Throat(6e-6,0), 7:Throat(7e-6,0), 8:Throat(8e-6,0), 9:Throat(9e-6,0)}, \
              3:{1:Throat(1e-6,3), 2:Throat(2e-6,3), 3:Throat(3e-6,3), 4:Throat(4e-6,3), 5:Throat(5e-6,3), 6:Throat(6e-6,3), 7:Throat(7e-6,3), 8:Throat(8e-6,3), 9:Throat(9e-6,3)}, \
              4:{1:Throat(1e-6,4), 2:Throat(2e-6,4), 3:Throat(3e-6,4), 4:Throat(4e-6,4), 5:Throat(5e-6,4), 6:Throat(6e-6,4), 7:Throat(7e-6,4), 8:Throat(8e-6,4), 9:Throat(9e-6,4)}, \
              5:{1:Throat(1e-6,5), 2:Throat(2e-6,5), 3:Throat(3e-6,5), 4:Throat(4e-6,5), 5:Throat(5e-6,5), 6:Throat(6e-6,5), 7:Throat(7e-6,5), 8:Throat(8e-6,5), 9:Throat(9e-6,5)}, \
              6:{1:Throat(1e-6,6), 2:Throat(2e-6,6), 3:Throat(3e-6,6), 4:Throat(4e-6,6), 5:Throat(5e-6,6), 6:Throat(6e-6,6), 7:Throat(7e-6,6), 8:Throat(8e-6,6), 9:Throat(9e-6,6)}}
TT=[15, 25, 35, 45, 55, 65, 75, 85, 95] # --------------------------------------------- Throat Type
DToGray={1e-6:225, -1e-6:220, 2e-6:200, -2e-6:195, 3e-6:175, -3e-6:170, 4e-6:150, -4e-6:145, 5e-6:125, -5e-6:120, \
         6e-6:100, -6e-6: 95, 7e-6: 75, -7e-6: 70, 8e-6: 50, -8e-6: 45, 9e-6: 25, -9e-6: 20, 0:0}
TToRGB ={}


print(ThroatChoice[5][9].D)
haha=[]
for tt in TT:
	haha.append(tt)
	haha.append(-tt)
Dia=np.zeros(len(haha), float)
for ha in haha:
	flip      =ha/abs(ha)
	row       =abs(ha)%10
	col       =abs(ha)//10
	Dia=ThroatChoice[row][col].D*flip
	print(ha, ThroatChoice[row][col].D*flip)

