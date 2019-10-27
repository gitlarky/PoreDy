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
TT=[[td, 5] for td in TD] # ------------------------------------------------------------Throat Type

#============================ Basic Functions =====================================================
# Pick one from a list ----------------------------------------------------------------------------
def Pick(List=[], Index=0, Method='Random'):
	if Method=='Random':
		return List[random.randint(0, len(List)-1)]
	elif Method=='Rotate':
		return List[Index%len(List)]

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
	def UseRandom(self, SubC=[]): # SubC: Sub Throat Index Collection
		Choice=self.ExistType()
		if SubC!=[]:
			Choice=set(SubC) & set(Choice)
		Index=Pick(List=Choice)
		self.ThroatCount[Index]-=1
		return self.ThroatType[Index]

	# Return one back to ThroatPool ---------------------------------------------------------------
	def Return(self, Index=0, T=[], Count=1, Method='ByIndex'):
		if Method=='ByIndex':
			Pick(List=self.ThroatCount, Index=Index, Method='Rotate')+=Count
			T=Pick(List=self.ThroatType, Index=Index, Method='Rotate')
		elif Method=='ByThroatType':
			for i in range(len(self.ThroatType)):
				if self.ThroatType[i]==T:
					Index=i
					self.ThroatCount[i]+=Count
		return True

#============================ Class PoreNetwork ===================================================
class PoreNetwork(object):
	# Initialize ----------------------------------------------------------------------------------
	def __init__(self, Name='PoreNetwork', Nx=20, Ny=20, Open=['N'], Cross=False, FixV=False, 
		         StrTType=[], CrsTType=[], StrTPool=ThroatPool(), CrsTPool=ThroatPool()):
		self.Name          =Name
		self.Nx            =Nx
		self.Ny            =Ny
		self.Open          =Open
		self.Cross         =Cross
		self.FixV          =FixV
		self.Mx            =2*self.Nx+1
		self.My            =2*self.Ny+1
		self.Matrix        =[[[0, 0] for i in range(self.Mx)] for j in range(self.My)]
		self.StrTCount     =self.Ny*(self.Nx+1)+self.Nx*(self.Ny+1)
		self.CrsTCount     =self.Nx*self.Ny if Cross else 0
		self.VTRange       =[[0, self.Nx], [1, self.Ny]]
		self.HTRange       =[[1, self.Nx], [0, self.Ny]]
		self.CTRange       =[[1, self.Nx], [1, self.Ny]]
		if 'W' in self.Open:
			self.VTRange[0][0]=1
			self.StrTCount   -=self.Ny
		if 'E' in self.Open:
			self.VTRange[0][1]=self.Nx-1
			self.StrTCount   -=self.Ny
		if 'S' in self.Open:
			self.HTRange[0][0]=1
			self.StrTCount   -=self.Nx
		if 'N' in self.Open:
			self.HTRange[0][1]=self.Ny-1
			self.StrTCount   -=self.Nx
		
		if (StrTType==[]) and (StrTPool.TotalCount!=0):
			self.StrTPool     =StrTPool
			self.StrTType     =StrTPool.ThroatType
		elif (StrTType!=[] and StrTPool.TotalCount==0):
			self.StrTType     =StrTType
			if self.FixV:
				AvgStrTN              =self.StrTCount//len(StrTType)
				RemStrTN              =self.StrTCount%len(StrTType)
				STN                   =[AvgStrTN for i in range(len(StrTType))]
				STN[len(StrTType)//2]+=RemStrTN
			else:
				STN           =[self.StrTCount for i in range(len(StrTType))]
			self.StrTPool     =ThroatPool(Name=self.Name+'.StrTPool', 
					                      ThroatType=self.StrTType, ThroatCount=STN)

		if (CrsTType==[]) and (CrsTPool.TotalCount!=0):
			self.CrsTPool     =CrsTPool
			self.CrsTType     =CrsTPool.ThroatType
		elif (CrsTType!=[] and CrsTPool.TotalCount==0):
			self.CrsTType     =CrsTType
			if self.FixV:
				AvgCrsTN              =self.CrsTCount//len(CrsTType)
				RemCrsTN              =self.CrsTCount%len(CrsTType)
				CTN                   =[AvgCrsTN for i in range(len(CrsTType))]
				CTN[len(CrsTType)//2]+=RemCrsTN
			else:
				CTN           =[self.CrsTCount for i in range(len(CrsTType))]
			self.CrsTPool     =ThroatPool(Name=self.Name+'.CrsTPool', 
				                          ThroatType=self.CrsTType, ThroatCount=CTN)

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
					theThroat=self.CrsTPool.UseRandom()
					self.Matrix[i][j][0]=theThroat[0]*Pick(List=[1, -1])
					self.Matrix[i][j][1]=theThroat[1]
				else:
					self.Matrix[i][j]=self.StrTPool.UseRandom()

		self.TotTCount=self.StrTCount+self.CrsTCount

	# set a new name ------------------------------------------------------------------------------
	def SetName(self, Name):
		self.Name=Name
		return True

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
	
	# Assign Throat at a certain coordinate in Matrix ---------------------------------------------
	def AssignMC(self, i, j, T): # i, j is the index in Matrix, i in range [0, Mx), j in range [0, My)
		self.Matrix[i][j][0]=T[0]
		self.Matrix[i][j][1]=T[1]
		return True

	# Get Information for a certain coordinate in Matrix ------------------------------------------
	def GetMC(self, i, j):
		if i>=0 and i<self.Mx and j>=0 and j<self.My:
			if(i%2==0 and j%2!=0):
				I=i/2
				J=(j+1)/2
				if (I>=self.VTRange[0][0] and I<=self.VTRange[0][1] and 
					J>=self.VTRange[1][0] and J<=self.VTRange[1][1]):
					return ['VT'  , I, J, self.Matrix[i][j]]
				else:
					return ['VT'  , I, J]
			elif(i%2!=0 and j%2==0):
				I=(i+1)/2
				J=j/2
				if (I>=self.HTRange[0][0] and I<=self.HTRange[0][1] and 
					J>=self.HTRange[1][0] and J<=self.HTRange[1][1]):
					return ['HT'  , I, J, self.Matrix[i][j]]
				else:
					return ['HT'  , I, J]
			elif(i%2!=0 and j%2!=0):
				I=(i+1)/2
				J=(j+1)/2
				if self.Cross:
					return ['CT'  , I, J, self.Matrix[i][j]]
				else:
					return ['CT'  , I, J]
			elif(i%2==0 and j%2==0):
				I=i/2
				J=j/2
				return ['Pore', I, J]
		else:
			return []

	# Assign Throat at a certain position ---------------------------------------------------------
	def AssignT(self, TP, I, J, T):
		if TP=='VT':
			self.AssignMC(2*I  , 2*J-1, T)
		elif TP=='HT':
			self.AssignMC(2*I-1, 2*J  , T)
		elif TP=='CT':
			self.AssignMC(2*I-1, 2*J-1, T)
		else:
			return False
		return True
	
	# Assign Throat at a certain position ---------------------------------------------------------
	def GetT(self, TP, I, J):
		if TP=='Pore':
			i=2*I
			j=2*J
			return [i, j]
		elif TP=='VT':
			i=2*I
			j=2*J-1
			if (I>=self.VTRange[0][0] and I<=self.VTRange[0][1] and 
				J>=self.VTRange[1][0] and J<=self.VTRange[1][1]):
				return [i, j, self.Matrix[i][j]]
			else:
				return [i, j]
		elif TP=='HT':
			i=2*I-1
			j=2*J
			if (I>=self.HTRange[0][0] and I<=self.HTRange[0][1] and 
				J>=self.HTRange[1][0] and J<=self.HTRange[1][1]):
				return [i, j, self.Matrix[i][j]]
			else:
				return [i, j]
		elif TP=='CT':
			i=2*I-1
			j=2*J-1
			return [i, j, self.Matrix[i][i]]
		else:
			return []
	
	# Return Throat at a certain position back to appropriate pool --------------------------------
	def ReturnT(self, TP, I, J):
		T=GetT(TP, I, J)
		if ((TP=='VT' or TP=='HT') and len(T)==3):
			self.StrTPool.Return(T=T[2])
			self.AssignT(TP, I, J, [0, 0])
		elif TP=='CT' and len(T)==3:
			self.CrsTPool.Return(T=[abs(T[2][0]), T[2][1]])
			self.AssignT(TP, I, J, [0, 0])
		else:
			return False
		return True

	# Restore all throats back to appropriate pool ------------------------------------------------
	def Restore(self, TPs=[]):
		if 'VT' in TPs:
			for I in range(self.VTRange[0][0], self.VTRange[0][1]+1, 1):
				for J in range(self.VTRange[1][0], self.VTRange[1][1]+1, 1):
					self.ReturnT('VT', I, J)
		if 'HT' in TPs:
			for I in range(self.HTRange[0][0], self.HTRange[0][1]+1, 1):
				for J in range(self.HTRange[1][0], self.HTRange[1][1]+1, 1):
					self.ReturnT('HT', I, J)
		if 'CT' in TPs:
			for I in range(self.CTRange[0][0], self.CTRange[0][1]+1, 1):
				for J in range(self.CTRange[1][0], self.CTRange[1][1]+1, 1):
					self.ReturnT('CT', I, J)
		return True

	# Assign Throat Box at a certain position -----------------------------------------------------
	def AssignBox(self, TP='',
		                Start=[0, 0], End=[0, 0], Band=[0, 0],
		                SubC=[], 
		                GradIJ=[0, 0], RepeIJ=[0, 0], JumpIJ=[0, 0],
		                GradB=0, RepeB=0, JumpB=0):
		if TP=='VT':
			Range=self.VTRange
		elif TP=='HT':
			Range=self.HTRange
		elif TP=='CT':
			Range=self.CTRange
		if Start[0]==0 and End[0]==0:
			Start[0]=Range[0][0]
			End  [0]=Range[0][1]
		if Start[1]==0 and End[1]==0:
			Start[1]=Range[1][0]+1
			End  [1]=Range[1][1]+1
		Count=0
		for I in range(Start[0], End[0], 1):
			for J in range(Start[1], End[1], 1):
				if not (Band[0]>0 and Band[0]>0 # if not (Band defined and in the band)
					    and I>=Start[0]+Band[0] and I<=End[0]-1-Band[0] 
					    and J>=Start[1]+Band[1] and J<=End[1]-1-Band[1]):
					self.AssignT(TP, I, J, T)
					Count+=1
		return Count

	# Slice a Region from original Matrix ---------------------------------------------------------
	def AssignRegion(self, TType='',
		                   IRange=[], IJump=0, IRep=1):
		Count=0
		return Count
#============================ Create Pore-Network Samples =========================================
def CreatePoreNetworkSamples(Nx=20, Ny=20, Folder=''):
	# Determine What and How Many is each type of Throat ------------------------------------------
	FixVStrNet =PoreNetwork(Name='FixVStrNet', Nx=Nx, Ny=Ny, FixV=True, 
		                    StrTType=TT)
	FixVCrsNet =PoreNetwork(Name='FixVCrsNet', Nx=Nx, Ny=Ny, Cross=True, FixV=True, 
		                    StrTType=TT, CrsTType=TT)
	RandStrNet =PoreNetwork(Name='FixVStrNet', Nx=Nx, Ny=Ny, 
		                    StrTType=TT)
	RandCrsNet =PoreNetwork(Name='FixVCrsNet', Nx=Nx, Ny=Ny, Cross=True, 
		                    StrTType=TT, CrsTType=TT)

	print('Straight Throat Number:', FixVCrsNet.StrTCount)
	print('Cross    Throat Number:', FixVCrsNet.CrsTCount)
	print('Total    Throat Number:', FixVCrsNet.TotTCount)

	# print('Straight Throat Pool:\n\t', FixVStrPool.ThroatType, '\n\t', FixVStrPool.ThroatCount)
	# print('Cross    Throat Pool:\n\t', FixVCrsPool.ThroatType, '\n\t', FixVCrsPool.ThroatCount)
	FixVStrNet.Write()
	FixVCrsNet.Write()
	RandStrNet.Write()
	RandCrsNet.Write()
	
	SampleIndex=0

	# Create some Pore-Network with only random straight Throats ----------------------------------
	# for vir in range(Nx):
	# 	for vjr in range(1, Ny, 1):
	# 		for hir in range(1, Nx, 1):
	# 			for hjr in range(0, Ny-1, 1):

	# for vrs in range(TypeN):
	# 	for vcs in range(TypeN):
	# 		for hrs in range(TypeN):
	# 			for hcs in range(TypeN):

	

	









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