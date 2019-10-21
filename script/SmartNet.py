#==================================================================================================
# Title      : SmartNet
# Description: Smart Predictiing and Designing of Pore-Network
# Author     : Zhenyu Xu; westlark@outlook.com
# Start Time : 2019.09.01
# License    : Apache 2.0
#==================================================================================================

#============================ Module Importation ==================================================
import random
# from __future__ import absolute_import, division, print_function, unicode_literals

# import tensorflow as tf
# from tensorflow import keras

# print(tf.__version__)

# import numpy as np

# import os
#==================================================================================================
#============================ Global Definition ===================================================
TD=[1e-6, 2e-6, 3e-6, 4e-6, 5e-6, 6e-6, 7e-6, 8e-6, 9e-6] # ------Throat Parameter is Diameter here
TS=[0, 3, 4, 5, 6] # ---------------------------------------------Throat Crosssection Polygon Shape

#============================ Basic Functions =====================================================
# Pick one from a list ----------------------------------------------------------------------------
def Pick(List=[], Index='-1'):
	if Index<0 or Index>=len(List):
		return List[random.randint(0, len(List)-1)]
	else:
		return List[Index]

# Create file name --------------------------------------------------------------------------------
def Name(Prefix='Sample', Digit=6, Index=0):
	si=str(Index)
	if Digit<len(si):
		Digit=len(si)
	for i in range(Digit-len(si)):
		si='0'+si
	return Prefix+si

#============================ Class ThroatPool ====================================================
class ThroatPool(object):
	# Initialize ----------------------------------------------------------------------------------
	def __init__(self, Name='ThroatPool', ThroatType=[], ThroatCount=[]):
		self.Name       =Name
		self.ThroatType =ThroatType
		self.ThroatCount=ThroatCount
		self.TypeCount  =len(self.ThroatType)
		self.TotalCount =0
		for count in self.ThroatCount:
			self.TotalCount+=count

	# Throat type count ---------------------------------------------------------------------------
	def TypeCount(self):
		return len(self.ThroatType)

	# Find out the index of existing throat
	def ExistType(self):
		EI=[]
		for i in range(len(self.ThroatCount)):
			if self.ThroatCount[i]>0:
				EI.append(i)
		return EI

	# Total throat available ----------------------------------------------------------------------
	def ExistCount(self):
		count=0
		for type in self.ExistType():
			count+=self.ThroatCount[type]
		return count

	# Use one from ThroatPool ---------------------------------------------------------------------
	def Use(self, Index, Count=1):
		if (self.ThroatCount[Index]-Count)>=0:
			self.ThroatCount[Index]-=Count
		return self.ThroatType[Index]

	# Use one randomly from ThroatPool ------------------------------------------------------------
	def UseRandom(self):
		Index=self.ExistType()[random.randint(0, len(self.ExistType()))]
		self.ThroatCount[Index]-=1
		return self.ThroatType[Index]


#============================ Class PoreNetwork ===================================================
class PoreNetwork(object):
	# Initialize ----------------------------------------------------------------------------------
	def __init__(self, Name='PoreNetwork', Nx=20, Ny=20, Open=['N'], 
		         DiameterChoice=[], ShapeChoice=[], Cross=False):
		self.Name          =Name
		self.Nx            =Nx
		self.Ny            =Ny
		self.Open          =Open
		self.DiameterChoice=DiameterChoice
		self.ShapeChoice   =ShapeChoice
		self.Mx            =2*self.Nx+1
		self.My            =2*self.Ny+1
		self.Matrix        =[[[0, 0] for i in range(self.Mx)] for j in range(self.My)]
		self.StraightTCount=0
		self.CrossTCount   =0
		self.VTRange       =[[0, self.Nx], [1, self.Ny]]
		self.HTRange       =[[1, self.Nx], [0, self.Ny]]
		self.CTRange       =[[1, self.Nx], [1, self.Ny]]
		if 'W' in self.Open:
			self.VTRange[0][0]=1
		if 'E' in self.Open:
			self.VTRange[0][1]=self.Nx-1
		if 'S' in self.Open:
			self.HTRange[0][0]=1
		if 'N' in self.Open:
			self.HTRange[0][1]=self.Ny-1
		for i in range(self.Mx):
			for j in range(self.My):
				if  (i%2==0 and j%2==0):
					continue
				elif('W' in self.Open and i==0):
					continue
				elif('E' in self.Open and i==(self.Mx-1)):
					continue
				elif('S' in self.Open and j==0):
					continue
				elif('N' in self.Open and j==(self.My-1)):
					continue
				elif(Cross==False and (i%2!=0 and j%2!=0)):
					continue
				elif(Cross==True and (i%2!=0 and j%2!=0)):
					self.Matrix[i][j][0]=Pick(self.DiameterChoice, -1)*Pick([1, -1], -1)
					self.Matrix[i][j][1]=Pick(self.ShapeChoice, -1)
					self.CrossTCount+=1
				else:
					self.Matrix[i][j][0]=Pick(self.DiameterChoice, -1)
					self.Matrix[i][j][1]=Pick(self.ShapeChoice, -1)
					self.StraightTCount+=1
		self.TotalTCount   =self.StraightTCount+self.CrossTCount

	# set a new name ------------------------------------------------------------------------------
	def SetName(self, Name):
		self.Name=Name
		return True

	# Assign Throat at a certain coordinate in Matrix ---------------------------------------------
	def AssignT(self, i, j, T): # i, j is the index in Matrix, i in range [0, Mx), j in range [0, My)
		self.Matrix[i][j][0]=T[0]
		self.Matrix[i][j][1]=T[1]
		return True

	# Assign Vertical Throat at a certain position ------------------------------------------------
	def AssignVT(self, I, J, T): # I, J is the index in Nx and Ny, I in range [0, Nx], J in range [1, Ny]
		self.AssignT(2*I, 2*J-1, T)
		return True
	# Assign Horizontal Throat at a certain position ----------------------------------------------
	def AssignHT(self, I, J, T): # I, J is the index in Nx and Ny, I in range [1, Nx], J in range [0, Ny]
		self.AssignT(2*I-1, 2*J, T)
		return True
	# Assign Cross Throat at a certain position ---------------------------------------------------
	def AssignCT(self, I, J, T): # I, J is the index in Nx and Ny, I in range [1, Nx], J in range [1, Ny]
		self.AssignT(2*I-1, 2*J-1, T)
		return True

	# Assign Vertical Throat Box at a certain position --------------------------------------------
	def AssignVTBox(self, T, IRange=[], JRange=[], IBand=0, JBand=0):
		if IRange==[]: 
			IRange=self.VTRange[0]
		if JRange==[]:
			JRange=self.VTRange[1]
		Count=0
		for I in range(IRange[0], IRange[1], 1):
			for J in range(JRange[0], JRange[1], 1):
				if (IBand>0 and JBand>0 
					and I>=IRange[0]+IBand and I<=IRange[1]-IBand 
					and J>=JRange[0]+JBand and J<=JRange[1]-JBand):
					continue
				else:
					self.AssignVT(I, J, T)
					Count+=1
		return Count
	# Assign Horizontal Throat Box at a certain position ------------------------------------------
	def AssignHTBox(self, T, IRange=[], JRange=[], IBand=0, JBand=0):
		if IRange==[]:
			IRange=self.HTRange[0]
		if JRange==[]:
			JRange=self.HTRange[1]
		Count=0
		for I in range(IRange[0], IRange[1], 1):
			for J in range(JRange[0], JRange[1], 1):
				if (IBand>0 and JBand>0 
					and I>=IRange[0]+IBand and I<=IRange[1]-IBand 
					and J>=JRange[0]+JBand and J<=JRange[1]-JBand):
					continue
				else:
					self.AssignHT(I, J, T)
					Count+=1
		return Count
	# Assign Cross Throat Box at a certain position -----------------------------------------------
	def AssignCTBox(self, T, IRange=[], JRange=[], IBand=0, JBand=0):
		if IRange==[]:
			IRange=self.CTRange[0]
		if JRange==[]:
			JRange=self.CTRange[1]
		Count=0
		for I in range(IRange[0], IRange[1], 1):
			for J in range(JRange[0], JRange[1], 1):
				if (IBand>0 and JBand>0 
					and I>=IRange[0]+IBand and I<=IRange[1]-IBand 
					and J>=JRange[0]+JBand and J<=JRange[1]-JBand):
					continue
				else:
					self.AssignCT(I, J, T)
					Count+=1
		return Count

	# Slice a Region from original Matrix ---------------------------------------------------------
	def AssignVTRegion(self, IRange=[], IJump=0, IRep=1)
	# Write PoreNetwork Data File -----------------------------------------------------------------
	def Write(self):
		with open(self.Name+'.at', 'w') as wat:
			for j in range(self.My-1, -1, -1):
				for i in range(self.Mx):
					wat.write('% 9.6e\t% 9d\t' % (self.Matrix[i][j][0], self.Matrix[i][j][1]))
				wat.write('\n')
		wat.close()
		return True

	# Write PoreNetwork Pixel File ----------------------------------------------------------------
	def Pixel(self):
		return True

#============================ Create Pore-Network Samples =========================================
def CreatePoreNetworkSamples(Nx=20, Ny=20, Folder=''):
	# Determine What and How Many is each type of Throat ------------------------------------------
	TypeN      =len(TD)
	FixVStrNet =PoreNetwork(Name='FixVStrNet', Nx=Nx, Ny=Ny, DiameterChoice=TD, ShapeChoice=[5])
	FixVCrsNet =PoreNetwork(Name='FixVCrsNet', Nx=Nx, Ny=Ny, DiameterChoice=TD, ShapeChoice=[5], 
		                    Cross=True)
	RandStrNet =PoreNetwork(Name='RandStrNet', Nx=Nx, Ny=Ny, DiameterChoice=TD, ShapeChoice=[5])
	RandCrsNet =PoreNetwork(Name='RandCrsNet', Nx=Nx, Ny=Ny, DiameterChoice=TD, ShapeChoice=[5], 
		                    Cross=True)
	StrTN      =FixVCrsNet.StraightTCount
	CrsTN      =FixVCrsNet.CrossTCount
	TotTN      =StrTN+CrsTN
	print('Straight Throat Number:', StrTN)
	print('Cross    Throat Number:', CrsTN)
	print('Total    Throat Number:', TotTN)
	AvgStrTN   =StrTN//TypeN
	RemStrTN   =StrTN%TypeN
	AvgCrsTN   =CrsTN//TypeN
	RemCrsTN   =CrsTN%TypeN
	STN        =[AvgStrTN for i in range(TypeN)]
	STN[4]    +=RemStrTN
	CTN        =[AvgCrsTN for i in range(TypeN)]
	CTN[4]    +=RemCrsTN
	FixVStrPool=ThroatPool(Name='FixVStrPool', ThroatType=[[td, 5] for td in TD], ThroatCount=STN)
	FixVCrsPool=ThroatPool(Name='FixVCrsPool', ThroatType=[[td, 5] for td in TD], ThroatCount=CTN)
	print('Straight Throat Pool:\n\t', FixVStrPool.ThroatType, '\n\t', FixVStrPool.ThroatCount)
	print('Cross    Throat Pool:\n\t', FixVCrsPool.ThroatType, '\n\t', FixVCrsPool.ThroatCount)
	# FixVStrNet.Write()
	# FixVCrsNet.Write()
	# RandStrNet.Write()
	# RandCrsNet.Write()
	
	SampleIndex=0

	# Create some Pore-Network with only random straight Throats ----------------------------------
	for vir in range(Nx):
		for vjr in range(1, Ny, 1):
			for hir in range(1, Nx, 1):
				for hjr in range(0, Ny-1, 1):

	for vrs in range(TypeN):
		for vcs in range(TypeN):
			for hrs in range(TypeN):
				for hcs in range(TypeN):

	

	









# #============================ Basic Functions =====================================================
# def create_model():
# 	model = tf.keras.models.Sequential([
# 		keras.layers.Dense(512, activation=tf.keras.activations.relu, input_shape=(784,)),
# 		keras.layers.Dropout(0.2),
# 		keras.layers.Dense(10, activation=tf.keras.activations.softmax)])

# 	model.compile(optimizer=tf.keras.optimizers.Adam(),
#                 loss=tf.keras.losses.sparse_categorical_crossentropy,
#                 metrics=['accuracy'])

# 	return model

# checkpoint_path = "training_1/cp.ckpt"
# checkpoint_dir = os.path.dirname(checkpoint_path)

# # Create checkpoint callback
# cp_callback = tf.keras.callbacks.ModelCheckpoint(checkpoint_path,
#                                                  save_weights_only=True,
#                                                  verbose=1)

# model = create_model()

# model.fit(train_images, train_labels,  epochs = 10,
#           validation_data = (test_images,test_labels),
#           callbacks = [cp_callback])  # pass callback to training

# # This may generate warnings related to saving the state of the optimizer.
# # These warnings (and similar warnings throughout this notebook)
# # are in place to discourage outdated usage, and can be ignored.


# model = create_model()

# loss, acc = model.evaluate(test_images, test_labels)
# print("Untrained model, accuracy: {:5.2f}%".format(100*acc))

# model.load_weights(checkpoint_path)
# loss,acc = model.evaluate(test_images, test_labels)
# print("Restored model, accuracy: {:5.2f}%".format(100*acc))


#============================ Main Program ========================================================
# CreatePoreNetworkSamples(Nx=20, Ny=20, Folder='/home/xu/work/PoreNetwork2020Samples')
# CreatePoreNetworkSamples(Nx=10, Ny=10, Folder='/home/xu/work/PoreNetwork1010Samples')
# CreatePoreNetworkSamples(Nx=40, Ny=40, Folder='/home/xu/work/PoreNetwork4040Samples')
CreatePoreNetworkSamples(Nx=3, Ny=3, Folder='/home/xu/work/PoreNetwork1010Samples')