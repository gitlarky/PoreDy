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
	def TypeCount():
		return len(self.ThroatType)

	# Find out the index of existing throat
	def ExistType():
		EI=[]
		for i in range(len(self.ThroatCount)):
			if self.ThroatCount[i]>0:
				EI.append(i)
		return EI

	# Total throat available ----------------------------------------------------------------------
	def ExistCount():
		count=0
		for type in ExistType():
			count+=self.ThroatCount[type]
		return count

	# Use one from ThroatPool ---------------------------------------------------------------------
	def Use(Index, Count=1):
		if (self.ThroatCount[Index]-Count)>=0:
			self.ThroatCount[Index]-=Count
		return self.ThroatType[Index]

#============================ Class PoreNetwork ===================================================
class PoreNetwork(object):
	# Initialize ----------------------------------------------------------------------------------
	def __init__(self, Name='PoreNetwork', Nx=20, Ny=20, Open=['N'], Cross=False,
		         DiameterChoice=TD, ShapeChoice=TS):
		self.Name          =Name
		self.Nx            =Nx
		self.Ny            =Ny
		self.Open          =Open
		self.DiameterChoice=DiameterChoice
		self.ShapeChoice   =ShapeChoice
		self.Mx            =2*self.Nx+1
		self.My            =2*self.Ny+1
		self.Matrix        =[[[0, 0] for i in range(Mx)] for j in range(My)]
		for i in range(self.Mx):
			for j in range(self.My):
				if  (i%2==0 and j%2==0):
				elif('W' in self.Open and i==0):
				elif('E' in self.Open and i==(self.Mx-1)):
				elif('S' in self.Open and j==0):
				elif('N' in self.Open and j==(self.My-1)):
				elif(Cross==False and (i%2!=0 and j%2!=0)):
				elif(Cross==True and (i%2!=0 and j%2!=0)):
					Matrix[i][j][0]=Pick(self.DiameterChoice, -1)*Pick([1, -1], -1)
					Matrix[i][j][1]=Pick(self.ShapeChoice, -1)
				else:
					Matrix[i][j][0]=Pick(self.DiameterChoice, -1)
					Matrix[i][j][1]=Pick(self.ShapeChoice, -1)

	# Assign Throat at a certain coordinate in Matrix ---------------------------------------------
	def AssignT(i, j, V): # ----i, j is the index in Matrix, i in range [0, Mx), j in range [0, My)
		self.Matrix[i][j][0]=V[0]
		self.Matrix[i][j][1]=V[1]
		return True

	# Assign Vertical Throat at a certain position ------------------------------------------------
	def AssignVT(I, J, V): # I, J is the index in Nx and Ny, I in range [0, Nx], J in range [1, Ny]
		AssignT(2*I, 2*J-1, V)
		return True
	# Assign Horizontal Throat at a certain position ----------------------------------------------
	def AssignHT(I, J, V): # I, J is the index in Nx and Ny, I in range [1, Nx], J in range [0, Ny]
		AssignT(2*I-1, 2*J, V)
		return True
	# Assign Cross Throat at a certain position ---------------------------------------------------
	def AssignCT(I, J, V): # I, J is the index in Nx and Ny, I in range [1, Nx], J in range [1, Ny]
		AssignT(2*I-1, 2*J-1, V)
		return True

	# Assign Vertical Throat Box at a certain position --------------------------------------------
	def AssignVTBox(V, IRange=[0, self.Nx], JRange=[1, self.Ny], Band=0):
		Count=0
		for I in range(IRange[0], IRange[1], 1):
			for J in range(JRange[0], JRange[1], 1):
				if (Band>0 and I>=IRange[0]+Band and I<=IRange[1]-Band 
					       and J>=JRange[0]+Band and J<=JRange[1]-Band):
				else:
					AssignVT(I, J, V)
					Count+=1
		return Count
	# Assign Horizontal Throat Box at a certain position ------------------------------------------
	def AssignHTBox(V, IRange=[1, self.Nx], JRange=[0, self.Ny], Band=0):
		Count=0
		for I in range(IRange[0], IRange[1], 1):
			for J in range(JRange[0], JRange[1], 1):
				if (Band>0 and I>=IRange[0]+Band and I<=IRange[1]-Band 
					       and J>=JRange[0]+Band and J<=JRange[1]-Band):
				else:
					AssignHT(I, J, V)
					Count+=1
		return Count
	# Assign Cross Throat Box at a certain position -----------------------------------------------
	def AssignCTBox(V, IRange=[1, self.Nx], JRange=[1, self.Ny], Band=0):
		Count=0
		for I in range(IRange[0], IRange[1], 1):
			for J in range(JRange[0], JRange[1], 1):
				if (Band>0 and I>=IRange[0]+Band and I<=IRange[1]-Band 
					       and J>=JRange[0]+Band and J<=JRange[1]-Band):
				else:
					AssignCT(I, J, V)
					Count+=1
		return Count

	# Write PoreNetwork Data File -----------------------------------------------------------------
	def Write():
		with open(self.Name+'.at', 'w') as wat:
			for i in range(self.Mx):
				for j in range(self.My):
					wat.write('% 9.6e\t% 9d\t' % (self.Matrix[i][j][0], self.Matrix[i][j][1]))
				wat.write('\n')
		wat.close()
		return True

	# Write PoreNetwork Pixel File ----------------------------------------------------------------
	def Pixel():
		return True

#============================ Create Pore-Network Samples =========================================

def CreatePoreNetworkSamples(Nx=0, Ny=0, Folder=''):
	# Determine What and How Many is each type of Throat ------------------------------------------
	TP         =[0, 1e-6, 2e-6, 3e-6, 4e-6, 5e-6, 6e-6, 7e-6, 8e-6, 9e-6] # Parameter is Diameter here
	PolyN      =[0, 3, 4, 5, 6]
	TypeN      =len(TP)
	print('Throat Parameter   :', TP)
	ShoTNum    =Ny*(Nx+1)+Ny*Nx # Short(Straight) Throats
	LonTNum    =Ny*Nx           # Long (Inclined) Throats
	TotTNum    =ShoTNum+LonTNum
	AvgShoTNum =ShoTNum//TypeN
	RemShoTNum =ShoTNum%TypeN
	STN        =[AvgShoTNum    for i in range(TypeN)]
	STN[4]    +=RemShoTNum
	print('Short Throat Number:', STN)
	AvgLonTNum =LonTNum//TypeN
	RemLonTNum =LonTNum%TypeN
	LTN        =[AvgLonTNum    for i in range(TypeN)]
	LTN[4]    +=RemLonTNum
	print('Long  Throat Number:', LTN)
	TN         =[STN[i]+LTN[i] for i in range(TypeN)]
	print('Total Throat Number:', TN)
	SampleIndex=0
	VRowSBunch =AvgShoTNum//(Nx+1)
	VColSBunch =AvgShoTNum//Ny
	HRowSBunch =AvgShoTNum//Nx
	HColSBunch =AvgShoTNum//Ny
	RowLBunch  =AvgLonTNum//Nx
	ColLBunch  =AvgLonTNum//Ny

	# Create some Pore-Network with only short Throats
	for offsetx in range(TypeN-1):
		for offsety in range(TypeN-1):
			for deltax in range(TypeN-1):
				for deltay in range(TypeN-1):
					for vr in range(VRowSBunch+1):
						for vc in range(VColSBunch+1):
							for hr in range(HRowSBunch+1):
								for hc in range(HColSBunch+1):
									Matrix=[[[0, 0] for i in range(Nx*2+1)] for j in range(Ny*2)]
									for j in range(1, 2*Ny, 2):
										addy+=deltay
										for i in range(0, 2*Nx+1, 2):
											Matrix[i][j][0]=TP(offsety+addy)
											Matrix[i][j][1]=PolyN

									for j in range(0, 2*Ny, 2):
										for i in range(1, 2*Nx, 2):
											Matrix[i][j][0]=TP(offsety+addy)
											Matrix[i][j][1]=PolyN
									
									WriteAT(Name(Prefix='Row', Index=SampleIndex), Matrix)
									SampleIndex+=1

	# Create some Pore-Network with both short and long Throats

	# for offset in range(Nx)

	# with open(CaseName+'.at', 'w') as at:
	# 	for j in range(Ny*2):
	# 		for i in range(Nx*2):
	# 			if(i%2==0 and j%2==0):
	# 				Diameter=0
	# 				PolyN   =0
	# 			elif():
	# 				Diameter=
	# 				PolyN   =5
	# 			at.write('% 9d\t% 9.6e\t' % (Diameter, PolyN))
	# 		at.write('\n')






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
CreatePoreNetworkSamples(Nx=20, Ny=20, Folder='/home/xu/work/PoreNetwork2020Samples')
# CreatePoreNetworkSamples(Nx=10, Ny=10, Folder='/home/xu/work/PoreNetwork1010Samples')
# CreatePoreNetworkSamples(Nx=40, Ny=40, Folder='/home/xu/work/PoreNetwork4040Samples')
# CreatePoreNetworkSamples(Nx=3, Ny=3, Folder='/home/xu/work/PoreNetwork1010Samples')