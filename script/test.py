
from pyvista import examples
import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
# X = np.random.random((100, 100)) # sample 2D array
# print(X)
# mpl.pyplot.imshow(X, cmap="gray")
# mpl.pyplot.show()
mesh = examples.download_st_helens().warp_by_scalar()
# Add scalar array with range (0, 100) taht correlates with elevation
mesh['values'] = pv.plotting.normalize(mesh['Elevation']) * 100
# Define the colors we want to use
blue = np.array([12/256, 238/256, 246/256, 1])
black = np.array([11/256, 11/256, 11/256, 1])
grey = np.array([189/256, 189/256, 189/256, 1])
yellow = np.array([255/256, 247/256, 0/256, 1])
red = np.array([1, 0, 0, 1])

mapping = np.linspace(mesh['values'].min(), mesh['values'].max(), 256)
newcolors = np.empty((256, 4))
newcolors[mapping >= 80] = red
newcolors[mapping < 80] = grey
newcolors[mapping < 55] = yellow
newcolors[mapping < 30] = blue
newcolors[mapping < 1] = black

# Make the colormap from the listed colors
my_colormap = ListedColormap(newcolors)



mesh.plot(scalars='values', cmap=my_colormap)
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