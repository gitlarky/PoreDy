
from pyvista import examples
import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import collections as co
Throat=co.namedtuple('T', 'D, S')
ThroatChoice={0:{1:Throat(1e-6,0), 2:Throat(2e-6,0), 3:Throat(3e-6,0), 4:Throat(4e-6,0), 5:Throat(5e-6,0), 6:Throat(6e-6,0), 7:Throat(7e-6,0), 8:Throat(8e-6,0), 9:Throat(9e-6,0)}\
             ,3:{1:Throat(1e-6,3), 2:Throat(2e-6,3), 3:Throat(3e-6,3), 4:Throat(4e-6,3), 5:Throat(5e-6,3), 6:Throat(6e-6,3), 7:Throat(7e-6,3), 8:Throat(8e-6,3), 9:Throat(9e-6,3)}\
             ,4:{1:Throat(1e-6,4), 2:Throat(2e-6,4), 3:Throat(3e-6,4), 4:Throat(4e-6,4), 5:Throat(5e-6,4), 6:Throat(6e-6,4), 7:Throat(7e-6,4), 8:Throat(8e-6,4), 9:Throat(9e-6,4)}\
             ,5:{1:Throat(1e-6,5), 2:Throat(2e-6,5), 3:Throat(3e-6,5), 4:Throat(4e-6,5), 5:Throat(5e-6,5), 6:Throat(6e-6,5), 7:Throat(7e-6,5), 8:Throat(8e-6,5), 9:Throat(9e-6,5)}\
             ,6:{1:Throat(1e-6,6), 2:Throat(2e-6,6), 3:Throat(3e-6,6), 4:Throat(4e-6,6), 5:Throat(5e-6,6), 6:Throat(6e-6,6), 7:Throat(7e-6,6), 8:Throat(8e-6,6), 9:Throat(9e-6,6)}}

ThroatChoice[5]
ThroatChoice[4][4]

# # List Operation ----------------------------------------------------------------------------------
# def SetUnion(Base=[], SubC=[]):
# 	final=[]
# 	for item in Base:
# 		final.append(item)
# 	for item in SubC:
# 		if not item in Base:
# 			final.append(item)
# 	return final

# # List Operation ----------------------------------------------------------------------------------
# def SetJoint(Base=[], SubC=[]):
# 	final=[]
# 	for item in Base:
# 		if item in SubC:
# 			final.append(item)
# 	return final

# # List Operation ----------------------------------------------------------------------------------
# def SetDifference(Base=[], ExC=[]):
# 	final=[]
# 	for item in Base:
# 		if item in ExC:
# 			continue
# 		else:
# 			final.append(item)
# 	return final

# def haha(t):
# 	a=t
# 	a=[0,0]
# 	return a



# a=[0,1,2,3,4,5,6,7,8,9]
# b=[0,2,4,6,8]
# c=[1,3,5,7,9]
# d=[3,4,5,6]
# print('1', SetUnion(a,b))
# print('2', SetUnion(b,c))
# print('3', SetUnion(c,d))
# print(a,b,c,d)
# # a=[0,1,2,3,4,5,6,7,8,9]
# # b=[0,2,4,6,8]
# # c=[1,3,5,7,9]
# # d=[3,4,5,6]
# print('4', SetJoint(a,b))
# print('5', SetJoint(b,c))
# print('6', SetJoint(c,d))
# print(a,b,c,d)
# # a=[0,1,2,3,4,5,6,7,8,9]
# # b=[0,2,4,6,8]
# # c=[1,3,5,7,9]
# # d=[3,4,5,6]
# print('7', SetDifference(a,b))
# print('8', SetDifference(b,c))
# print('9', SetDifference(c,d))
# print(a,b,c,d)

# print('10', haha(b))
# print('10', b)