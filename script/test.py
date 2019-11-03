import matplotlib.pyplot as plt
import numpy as np
 
# X = np.random.random((100, 100)) # sample 2D array
# print(X)
# mpl.pyplot.imshow(X, cmap="gray")
# mpl.pyplot.show()


# List Operation ----------------------------------------------------------------------------------
def SetUnion(Base=[], SubC=[]):
	final=[]
	for item in Base:
		final.append(item)
	for item in SubC:
		if not item in Base:
			final.append(item)
	return final

# List Operation ----------------------------------------------------------------------------------
def SetJoint(Base=[], SubC=[]):
	final=[]
	for item in Base:
		if item in SubC:
			final.append(item)
	return final

# List Operation ----------------------------------------------------------------------------------
def SetDifference(Base=[], ExC=[]):
	final=[]
	for item in Base:
		if item in ExC:
			continue
		else:
			final.append(item)
	return final





a=[0,1,2,3,4,5,6,7,8,9]
b=[0,2,4,6,8]
c=[1,3,5,7,9]
d=[3,4,5,6]
print('1', SetUnion(a,b))
print('2', SetUnion(b,c))
print('3', SetUnion(c,d))
print(a,b,c,d)
# a=[0,1,2,3,4,5,6,7,8,9]
# b=[0,2,4,6,8]
# c=[1,3,5,7,9]
# d=[3,4,5,6]
print('4', SetJoint(a,b))
print('5', SetJoint(b,c))
print('6', SetJoint(c,d))
print(a,b,c,d)
# a=[0,1,2,3,4,5,6,7,8,9]
# b=[0,2,4,6,8]
# c=[1,3,5,7,9]
# d=[3,4,5,6]
print('7', SetDifference(a,b))
print('8', SetDifference(b,c))
print('9', SetDifference(c,d))
print(a,b,c,d)