#==================================================================================================
# Title      : SmartNet
# Description: Smart Predictiing and Designing of Pore-Network
# Author     : Zhenyu Xu; westlark@outlook.com
# Start Time : 2019.09.01
# License    : Apache 2.0
#==================================================================================================

#============================ Module Importation ==================================================
import random
import pickle
import subprocess
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
def Name(Prefix='Sample', Index=0, Digit=6):
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
	def ExistTypeCount(self):
		etc=0
		for c in self.ThroatCount:
			if c>0:
				etc+=1
		return etc

	# Find out the index of existing throat
	def ExistTypeIndex(self):
		eti=[]
		for i in range(len(self.ThroatCount)):
			if self.ThroatCount[i]>0:
				eti.append(i)
		return eti

	# Total throat available ----------------------------------------------------------------------
	def ExistThroatCount(self):
		count=0
		for type in self.ExistTypeIndex():
			count+=self.ThroatCount[type]
		return count

	# Use appointed one from ThroatPool -----------------------------------------------------------
	def UseAppoint(self, Index, Count=1):
		if (self.ThroatCount[Index]-Count)>=0:
			self.ThroatCount[Index]-=Count
			return self.ThroatType[Index]
		else:
			print('Not enough throats of this type, nothing assigned!')
			return [0, 0]

	# Use one randomly from ThroatPool ------------------------------------------------------------
	def UseRandom(self, SubC=[]): # SubC: Sub Throat Index Collection
		Choice=[]
		if SubC==[]:
			Choice=self.ExistTypeIndex()
		else:
			for eti in self.ExistTypeIndex():
				if eti in SubC:
					Choice.append(eti)
		Index=Pick(List=Choice)
		self.ThroatCount[Index]-=1
		return self.ThroatType[Index]

	# Use one appointed or choose randomly from pool exclude a sub --------------------------------
	def Use(self, Index, ExC=[]):
		if self.ThroatCount[Index]>=1:
			self.ThroatCount[Index]-=1
		else:
			Available=[]
			for eti in self.ExistTypeIndex():
				if eti in ExC:
					continue
				else:
					Available.append(eti)
			Index=Pick(List=Available)
			self.ThroatCount[Index]-=1
		return self.ThroatType[Index]

	# Return one back to ThroatPool ---------------------------------------------------------------
	def Return(self, Index=0, T=[], Count=1, Method='ByIndex'):
		print(T)
		if Method=='ByIndex':
			I=Pick(List=self.ThroatCount, Index=Index, Method='Rotate')
			self.ThroatCount[I]+=Count
			T=self.ThroatType[I]
			return True
		elif Method=='ByThroatType':
			for i in range(len(self.ThroatType)):
				if self.ThroatType[i]==T:
					Index=i
					self.ThroatCount[Index]+=Count
					return True
			return False
		else:
			return False

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
		self.CrsTCount     =self.Nx*self.Ny if self.Cross else 0
		self.VTRange       =[[0, self.Nx], [1, self.Ny]]
		self.HTRange       =[[1, self.Nx], [0, self.Ny]]
		self.CTRange       =[[1, self.Nx], [1, self.Ny]] if self.Cross else []
		if 'W' in self.Open:
			self.VTRange[0][0]=1
			self.StrTCount   -=self.Ny
		if 'E' in self.Open:
			self.VTRange[0][1]=self.Nx-1
			self.StrTCount   -=self.Ny
		if 'S' in self.Open:
			self.HTRange[1][0]=1
			self.StrTCount   -=self.Nx
		if 'N' in self.Open:
			self.HTRange[1][1]=self.Ny-1
			self.StrTCount   -=self.Nx
		
		if (StrTType==[]) and (StrTPool.TotalCount!=0):
			self.StrTPool     =StrTPool
			self.StrTType     =StrTPool.ThroatType
		elif (StrTType!=[] and StrTPool.TotalCount==0):
			self.StrTType     =StrTType
			STN=[]
			if self.FixV:
				AvgStrTN              =self.StrTCount//len(StrTType)
				RemStrTN              =self.StrTCount%len(StrTType)
				STN                   =[AvgStrTN for i in range(len(StrTType))]
				STN[len(StrTType)//2]+=RemStrTN
				# print(AvgStrTN, RemStrTN, STN)
			else:
				STN           =[self.StrTCount for i in range(len(StrTType))]
			self.StrTPool     =ThroatPool(Name=self.Name+'.StrTPool', 
					                      ThroatType=self.StrTType, ThroatCount=STN)
		else:
			self.StrTPool=ThroatPool()
			self.StrTType=[]

		if (CrsTType==[]) and (CrsTPool.TotalCount!=0):
			self.CrsTPool     =CrsTPool
			self.CrsTType     =CrsTPool.ThroatType
		elif (CrsTType!=[] and CrsTPool.TotalCount==0):
			self.CrsTType     =CrsTType
			CTN=[]
			if self.FixV:
				AvgCrsTN              =self.CrsTCount//len(CrsTType)
				RemCrsTN              =self.CrsTCount%len(CrsTType)
				CTN                   =[AvgCrsTN for i in range(len(CrsTType))]
				CTN[len(CrsTType)//2]+=RemCrsTN
				# print(AvgCrsTN, RemCrsTN, CTN)
			else:
				CTN           =[self.CrsTCount for i in range(len(CrsTType))]
			self.CrsTPool     =ThroatPool(Name=self.Name+'.CrsTPool', 
				                          ThroatType=self.CrsTType, ThroatCount=CTN)
		else:
			self.CrsTPool=ThroatPool()
			self.CrsTType=[]

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
				elif(self.Cross==False and (i%2!=0 and j%2!=0)):
					continue
				elif(self.Cross==True  and (i%2!=0 and j%2!=0)):
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

	# Report information --------------------------------------------------------------------------
	def Report(self):
		print(self.Name, ':')
		print('\tCross=', self.Cross, 'FixV=', self.FixV, 'Open: ', self.Open)
		print('\tNx=', self.Nx, 'Ny=', self.Ny, 'Mx=', self.Mx, 'My=', self.My)
		print('\tStrTCount=', self.StrTCount, 'CrsTCount=', self.CrsTCount, 'TotTCount=', self.TotTCount)
		print('\tVTRange=', self.VTRange, 'HTRange=', self.HTRange, 'CTRange=', self.CTRange)
		print('\tStrTType=', self.StrTType, 'StrTPool: ', self.StrTPool.ThroatType, self.StrTPool.ThroatCount)
		print('\tCrsTType=', self.CrsTType, 'CrsTPool: ', self.CrsTPool.ThroatType, self.CrsTPool.ThroatCount)

		return True
	# Write PoreNetwork Data File -----------------------------------------------------------------
	def Write(self, Folder=''):
		with open(Folder+self.Name+'.at', 'w') as wat:
			for j in range(self.My-1, -1, -1):
				for i in range(self.Mx):
					wat.write('% 9.6e\t% 9d\t' % (self.Matrix[i][j][0], self.Matrix[i][j][1]))
				wat.write('\n')
		wat.close()
		return True

	# Write PoreNetwork Pixel File ----------------------------------------------------------------
	def Pixel(self, Scale=1):
		return True
	
	# Dump the foil into a binary file ------------------------------------------------------------
	def Dump(self, Destination='', Disp=False):
		if not Destination: Destination=self.Name+'.net'
		pickle.dump(self, open(Destination, 'wb'))
		if Disp: print('Successfully saved pore-network file.')

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
		if   TP=='VT' and I>=self.VTRange[0][0] and I<=self.VTRange[0][1] and J>=self.VTRange[1][0] and J<=self.VTRange[1][1]:
			self.AssignMC(2*I  , 2*J-1, T)
		elif TP=='HT' and I>=self.HTRange[0][0] and I<=self.HTRange[0][1] and J>=self.HTRange[1][0] and J<=self.HTRange[1][1]:
			self.AssignMC(2*I-1, 2*J  , T)
		elif TP=='CT' and I>=self.CTRange[0][0] and I<=self.CTRange[0][1] and J>=self.CTRange[1][0] and J<=self.CTRange[1][1] and self.Cross:
			self.AssignMC(2*I-1, 2*J-1, T)
		else:
			return False
		return True

	# Judge if it is assigned or not --------------------------------------------------------------
	def Assigned(self, TP, I, J):
		if   TP=='VT' and self.Matrix[2*I  ][2*J-1] in self.StrTPool.ThroatType:
			return True
		elif TP=='HT' and self.Matrix[2*I-1][2*J  ] in self.StrTPool.ThroatType:
			return True
		elif TP=='CT' and self.Matrix[2*I-1][2*J-1] in self.CrsTPool.ThroatType:
			return True
		else:
			return False

	# Get Throat info at a certain position -------------------------------------------------------
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
			if self.Cross and I>=self.CTRange[0][0] and I<=self.CTRange[0][1] and J>=self.CTRange[1][0] and J<=self.CTRange[1][1]:
				return [i, j, self.Matrix[i][j]]
			else:
				return [i, j]
		else:
			return []
	
	# Return Throat at a certain position back to appropriate pool --------------------------------
	def ReturnT(self, TP, I, J):
		T=self.GetT(TP, I, J)
		print(T[2])
		if ((TP=='VT' or TP=='HT') and len(T)==3):
			self.StrTPool.Return(T=T[2], Method='ByThroatType')
			self.AssignT(TP, I, J, [0, 0])
		elif TP=='CT' and len(T)==3:
			self.CrsTPool.Return(T=[abs(T[2][0]), T[2][1]], Method='ByThroatType')
			self.AssignT(TP, I, J, [0, 0])
		else:
			return False
		return True

	# Restore all throats back to appropriate pool ------------------------------------------------
	def Restore(self, TPs=['VT', 'HT', 'CT']):
		if 'VT' in TPs:
			for I in range(self.VTRange[0][0], self.VTRange[0][1]+1, 1):
				for J in range(self.VTRange[1][0], self.VTRange[1][1]+1, 1):
					self.ReturnT('VT', I, J)
		if 'HT' in TPs:
			for I in range(self.HTRange[0][0], self.HTRange[0][1]+1, 1):
				for J in range(self.HTRange[1][0], self.HTRange[1][1]+1, 1):
					self.ReturnT('HT', I, J)
		if 'CT' in TPs and self.Cross:
			for I in range(self.CTRange[0][0], self.CTRange[0][1]+1, 1):
				for J in range(self.CTRange[1][0], self.CTRange[1][1]+1, 1):
					self.ReturnT('CT', I, J)
		return True

	# Assign Throat distribution in a box of certain position -------------------------------------
	def AssignBox(self, TP='',
		                Start=[0, 0], End=[0, 0], Band=[0, 0],
		                SubC=[], Grad=[1, 1, 0], Repeat=[1, 1, 1], Flip=0):
		if   TP=='VT':
			Range=self.VTRange
			TypeChoice=self.StrTType
		elif TP=='HT':
			Range=self.HTRange
			TypeChoice=self.StrTType
		elif TP=='CT' and self.Cross:
			Range=self.CTRange
			TypeChoice=self.CrsTType
		if Start[0]==0 and End[0]==0:
			Start[0]=Range[0][0]
			End  [0]=Range[0][1]
		if Start[1]==0 and End[1]==0:
			Start[1]=Range[1][0]+1
			End  [1]=Range[1][1]+1
		if SubC==[]:
			SubC=[i for i in range(len(TypeChoice))]
			ExC =[]
		else:
			ExC =SubC

		Count=0
		for I in range(Start[0], End[0], 1):
			for J in range(Start[1], End[1], 1):
				if Band==[0, 0] or (
				   (Band[0]> 0 and Band[1]> 0) and not(I>=Start[0]+Band[0] and I<=End[0]-Band[0]-1 and 
				                                       J>=Start[1]+Band[1] and J<=End[1]-Band[1]-1)) or (
				   (Band[0]> 0 and Band[1]==0) and not(I>=Start[0]+Band[0] and I<=End[0]-Band[0]-1)) or (
				   (Band[0]==0 and Band[1]> 0) and not(J>=Start[1]+Band[1] and J<=End[1]-Band[1]-1)):
					if self.Assigned(TP, I, J):
						self.ReturnT(TP, I, J)

					for ring in range(min(Band[0], Band[1])+1):
						if not (I>=Start[0]+ring+1 and I<=End[0]-ring-2 and 
							J>=Start[1]+ring+1 and J<=End[1]-ring-2) and (
							I>=Start[0]+ring   and I<=End[0]-ring-1 and 
							J>=Start[1]+ring   and J<=End[1]-ring-1):
							RingIndex=ring
					Index=Pick(List=SubC, Index=(I-Start[0])//Repeat[0]*Grad[0]+
						                        (J-Start[1])//Repeat[1]*Grad[1]+
						                        RingIndex   //Repeat[2]*Grad[2], Method='Rotate')
					if TP=='VT' or TP=='HT':
						T=self.StrTPool.Use(Index=Index, ExC=ExC)
					elif TP=='CT':
						T=self.CrsTPool.Use(Index=Index, ExC=ExC)
						if Flip==0:
							T=[Pick(List=[-1, 1])*T[0], T[1]]
						elif Flip==2 and I%2==0:
							T=[      -1     *T[0], T[1]]
						elif Flip==3 and J%2==0:
							T=[      -1     *T[0], T[1]]
						else: # Flip= -1 or 1
							T=[Flip         *T[0], T[1]]
					self.AssignT(TP, I, J, T)
					Count+=1
		
		return Count

	# Slice a Region from original Matrix ---------------------------------------------------------
	def RandomRest(self):
		Count=0
		for I in range(self.VTRange[0][0], self.VTRange[0][1]+1, 1):
			for J in range(self.VTRange[1][0], self.VTRange[1][1]+1, 1):
				if not self.Assigned('VT', I, J):
					self.AssignT('VT', I, J, self.StrTPool.UseRandom())
					Count+=1
		for I in range(self.HTRange[0][0], self.HTRange[0][1]+1, 1):
			for J in range(self.HTRange[1][0], self.HTRange[1][1]+1, 1):
				if not self.Assigned('HT', I, J):
					self.AssignT('HT', I, J, self.StrTPool.UseRandom())
					Count+=1
		if self.Cross:
			for I in range(self.CTRange[0][0], self.CTRange[0][1]+1, 1):
				for J in range(self.CTRange[1][0], self.CTRange[1][1]+1, 1):
					if not self.Assigned('CT', I, J):
						self.AssignT('CT', I, J, self.CrsTPool.UseRandom())
						Count+=1

		return Count

#============================ Create Pore-Network Samples =========================================
def CreatePoreNetworkSamples(Nx=20, Ny=20, Folder=''):
	# Prepare the folders -------------------------------------------------------------------------
	cmd='mkdir '+Folder
	subprocess.call(cmd, shell=True)
	cmd='mkdir '+Folder+'RandStrNet/'
	subprocess.call(cmd, shell=True)
	cmd='mkdir '+Folder+'RandCrsNet/'
	subprocess.call(cmd, shell=True)
	cmd='mkdir '+Folder+'FixVStrNet/'
	subprocess.call(cmd, shell=True)
	cmd='mkdir '+Folder+'FixVCrsNet/'
	subprocess.call(cmd, shell=True)
	# Initialize 4 kinds of networks --------------------------------------------------------------
	RandStrNet =PoreNetwork(Name='RandStrNet', Nx=Nx, Ny=Ny, Open=['N'], Cross=False, FixV=False,
		                    StrTType=TT)
	RandCrsNet =PoreNetwork(Name='RandCrsNet', Nx=Nx, Ny=Ny, Open=['N'], Cross=True , FixV=False,
		                    StrTType=TT, CrsTType=TT)
	FixVStrNet =PoreNetwork(Name='FixVStrNet', Nx=Nx, Ny=Ny, Open=['N'], Cross=False, FixV=True , 
		                    StrTType=TT)
	FixVCrsNet =PoreNetwork(Name='FixVCrsNet', Nx=Nx, Ny=Ny, Open=['N'], Cross=True , FixV=True , 
		                    StrTType=TT, CrsTType=TT)
	RandStrNet.Report()
	RandCrsNet.Report()
	FixVStrNet.Report()
	FixVCrsNet.Report()

	RandStrNet.Restore()
	RandCrsNet.Restore()
	FixVStrNet.Restore()
	FixVCrsNet.Restore()

	RandStrNet.Report()
	RandCrsNet.Report()
	FixVStrNet.Report()
	FixVCrsNet.Report()


	# Generate Random Volume Network --------------------------------------------------------------
	SampleIndex=0
	for B in [0, 4, 8]:
		for SC in [[0, 1, 2], [4, 3, 5], [8, 7, 6], [2, 4, 6, 8], [7, 5, 3, 1]]:
			for IG in [0, 1, 2]:
				for JG in [0, 1, 2]:
					for IR in [1, 2, 3, 4, 5]:
						for JR in [1, 2, 3, 4, 5]:
							IS=0
							JS=0
							IE=0
							JE=0
							IB=B
							JB=B
							RG=0
							RR=1
							# ---------------------------------------------------------------------
							RandStrNet.Restore()
							if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<8:
								RandStrNet.AssignBox(TP='VT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
									                 SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
							RandStrNet.RandomRest()
							RandStrNet.SetName(Name=Name(Prefix='RandStrNet', Index=SampleIndex))
							RandStrNet.Write(Folder+'RandStrNet/')

							RandCrsNet.Restore()
							if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<8:
								RandCrsNet.AssignBox(TP='VT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
							CrsPick=Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
							if CrsPick<2:
								RandCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 2)
							elif CrsPick>7:
								RandCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 3)
							elif CrsPick==3 or CrsPick==4:
								RandCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip=-1)
							elif CrsPick==5 or CrsPick==6:
								RandCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 1)
							else:
								RandCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 0)
							RandCrsNet.RandomRest()
							RandCrsNet.SetName(Name=Name(Prefix='RandCrsNet', Index=SampleIndex))
							RandCrsNet.Write(Folder+'RandCrsNet/')
							SampleIndex+=1
							# ---------------------------------------------------------------------

							# ---------------------------------------------------------------------
							RandStrNet.Restore()
							if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<8:
								RandStrNet.AssignBox(TP='HT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
									                 SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
							RandStrNet.RandomRest()
							RandStrNet.SetName(Name=Name(Prefix='RandStrNet', Index=SampleIndex))
							RandStrNet.Write(Folder+'RandStrNet/')

							RandCrsNet.Restore()
							if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<8:
								RandCrsNet.AssignBox(TP='HT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
									                 SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
							CrsPick=Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
							if CrsPick<2:
								RandCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 2)
							elif CrsPick>7:
								RandCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 3)
							elif CrsPick==3 or CrsPick==4:
								RandCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip=-1)
							elif CrsPick==5 or CrsPick==6:
								RandCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 1)
							else:
								RandCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 0)
							RandCrsNet.RandomRest()
							RandCrsNet.SetName(Name=Name(Prefix='RandCrsNet', Index=SampleIndex))
							RandCrsNet.Write(Folder+'RandCrsNet/')
							SampleIndex+=1
							# ---------------------------------------------------------------------

							# ---------------------------------------------------------------------
							RandStrNet.Restore()
							if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<1:
								RandStrNet.AssignBox(TP='VT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
									                 SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
							if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<1:
								RandStrNet.AssignBox(TP='HT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
									                 SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
							RandStrNet.RandomRest()
							RandStrNet.SetName(Name=Name(Prefix='RandStrNet', Index=SampleIndex))
							RandStrNet.Write(Folder+'RandStrNet/')

							RandCrsNet.Restore()
							if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<1:
								RandCrsNet.AssignBox(TP='VT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
									                 SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
							if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<1:
								RandCrsNet.AssignBox(TP='HT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
									                 SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
							CrsPick=Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
							if CrsPick<2:
								RandCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 2)
							elif CrsPick>7:
								RandCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 3)
							elif CrsPick==3 or CrsPick==4:
								RandCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip=-1)
							elif CrsPick==5 or CrsPick==6:
								RandCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 1)
							else:
								RandCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 0)
							RandCrsNet.RandomRest()
							RandCrsNet.SetName(Name=Name(Prefix='RandCrsNet', Index=SampleIndex))
							RandCrsNet.Write(Folder+'RandCrsNet/')
							SampleIndex+=1
							# ---------------------------------------------------------------------
	print(SampleIndex)
	# ---------------------------------------------------------------------------------------------

	# Generate Fixed Volume Network ---------------------------------------------------------------
	SampleIndex=0
	for B in [0, 2, 4, 6, 8]:
		for SC in [[0, 1, 2, 3, 4], [8, 7, 6, 5, 4], [4, 0, 8, 2, 6], [4, 1, 7, 3, 5], [4, 5, 3, 6, 2]]:
			for IG in [1, 2, 3, 4]:
				for JG in [1, 2, 3, 4]:
					for IR in [1, 2, 3]:
						for JR in [1, 2, 3]:
							IS=0
							JS=0
							IE=0
							JE=0
							IB=B
							JB=B
							RG=0
							RR=1
							# ---------------------------------------------------------------------
							FixVStrNet.Restore()
							if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<8:
								FixVStrNet.AssignBox(TP='VT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
									                 SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
							FixVStrNet.RandomRest()
							FixVStrNet.SetName(Name=Name(Prefix='FixVStrNet', Index=SampleIndex))
							FixVStrNet.Write(Folder+'FixVStrNet/')

							FixVCrsNet.Restore()
							if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<8:
								FixVCrsNet.AssignBox(TP='VT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
							CrsPick=Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
							if CrsPick<2:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 2)
							elif CrsPick>7:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 3)
							elif CrsPick==3 or CrsPick==4:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip=-1)
							elif CrsPick==5 or CrsPick==6:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 1)
							else:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 0)
							FixVCrsNet.RandomRest()
							FixVCrsNet.SetName(Name=Name(Prefix='FixVCrsNet', Index=SampleIndex))
							FixVCrsNet.Write(Folder+'FixVCrsNet/')
							SampleIndex+=1
							# ---------------------------------------------------------------------

							# ---------------------------------------------------------------------
							FixVStrNet.Restore()
							if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<8:
								FixVStrNet.AssignBox(TP='HT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
									                 SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
							FixVStrNet.RandomRest()
							FixVStrNet.SetName(Name=Name(Prefix='FixVStrNet', Index=SampleIndex))
							FixVStrNet.Write(Folder+'FixVStrNet/')

							FixVCrsNet.Restore()
							if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<8:
								FixVCrsNet.AssignBox(TP='HT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
									                 SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
							CrsPick=Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
							if CrsPick<2:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 2)
							elif CrsPick>7:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 3)
							elif CrsPick==3 or CrsPick==4:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip=-1)
							elif CrsPick==5 or CrsPick==6:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 1)
							else:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 0)
							FixVCrsNet.RandomRest()
							FixVCrsNet.SetName(Name=Name(Prefix='FixVCrsNet', Index=SampleIndex))
							FixVCrsNet.Write(Folder+'FixVCrsNet/')
							SampleIndex+=1
							# ---------------------------------------------------------------------

							# ---------------------------------------------------------------------
							FixVStrNet.Restore()
							if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<1:
								FixVStrNet.AssignBox(TP='VT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
									                 SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
							if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<1:
								FixVStrNet.AssignBox(TP='HT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
									                 SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
							FixVStrNet.RandomRest()
							FixVStrNet.SetName(Name=Name(Prefix='FixVStrNet', Index=SampleIndex))
							FixVStrNet.Write(Folder+'FixVStrNet/')

							FixVCrsNet.Restore()
							if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<1:
								FixVCrsNet.AssignBox(TP='VT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
									                 SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
							if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<1:
								FixVCrsNet.AssignBox(TP='HT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
									                 SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
							CrsPick=Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
							if CrsPick<2:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 2)
							elif CrsPick>7:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 3)
							elif CrsPick==3 or CrsPick==4:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip=-1)
							elif CrsPick==5 or CrsPick==6:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 1)
							else:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 0)
							FixVCrsNet.RandomRest()
							FixVCrsNet.SetName(Name=Name(Prefix='FixVCrsNet', Index=SampleIndex))
							FixVCrsNet.Write(Folder+'FixVCrsNet/')
							SampleIndex+=1
							# ---------------------------------------------------------------------
	print(SampleIndex)
	for B in [1, 3, 5, 7, 9]:
		for SC in [[0, 1, 2, 3, 4], [8, 7, 6, 5, 4], [4, 0, 8, 2, 6], [4, 1, 7, 3, 5], [4, 5, 3, 6, 2]]:
			for RG in [1, 2, 3, 4]:
							IB=B
							JB=B
							IS=0
							JS=0
							IE=0
							JE=0
							IG=0
							JG=0
							IR=1
							JR=1
							RR=1
							# ---------------------------------------------------------------------
							FixVStrNet.Restore()
							if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<8:
								FixVStrNet.AssignBox(TP='VT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
									                 SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
							FixVStrNet.RandomRest()
							FixVStrNet.SetName(Name=Name(Prefix='FixVStrNet', Index=SampleIndex))
							FixVStrNet.Write(Folder+'FixVStrNet/')

							FixVCrsNet.Restore()
							if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<8:
								FixVCrsNet.AssignBox(TP='VT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
							CrsPick=Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
							if CrsPick<2:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 2)
							elif CrsPick>7:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 3)
							elif CrsPick==3 or CrsPick==4:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip=-1)
							elif CrsPick==5 or CrsPick==6:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 1)
							else:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 0)
							FixVCrsNet.RandomRest()
							FixVCrsNet.SetName(Name=Name(Prefix='FixVCrsNet', Index=SampleIndex))
							FixVCrsNet.Write(Folder+'FixVCrsNet/')
							SampleIndex+=1
							# ---------------------------------------------------------------------

							# ---------------------------------------------------------------------
							FixVStrNet.Restore()
							if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<8:
								FixVStrNet.AssignBox(TP='HT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
									                 SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
							FixVStrNet.RandomRest()
							FixVStrNet.SetName(Name=Name(Prefix='FixVStrNet', Index=SampleIndex))
							FixVStrNet.Write(Folder+'FixVStrNet/')

							FixVCrsNet.Restore()
							if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<8:
								FixVCrsNet.AssignBox(TP='HT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
									                 SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
							CrsPick=Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
							if CrsPick<2:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 2)
							elif CrsPick>7:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 3)
							elif CrsPick==3 or CrsPick==4:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip=-1)
							elif CrsPick==5 or CrsPick==6:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 1)
							else:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 0)
							FixVCrsNet.RandomRest()
							FixVCrsNet.SetName(Name=Name(Prefix='FixVCrsNet', Index=SampleIndex))
							FixVCrsNet.Write(Folder+'FixVCrsNet/')
							SampleIndex+=1
							# ---------------------------------------------------------------------

							# ---------------------------------------------------------------------
							FixVStrNet.Restore()
							if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<1:
								FixVStrNet.AssignBox(TP='VT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
									                 SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
							if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<1:
								FixVStrNet.AssignBox(TP='HT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
									                 SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
							FixVStrNet.RandomRest()
							FixVStrNet.SetName(Name=Name(Prefix='FixVStrNet', Index=SampleIndex))
							FixVStrNet.Write(Folder+'FixVStrNet/')

							FixVCrsNet.Restore()
							if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<1:
								FixVCrsNet.AssignBox(TP='VT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
									                 SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
							if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<1:
								FixVCrsNet.AssignBox(TP='HT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
									                 SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
							CrsPick=Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
							if CrsPick<2:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 2)
							elif CrsPick>7:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 3)
							elif CrsPick==3 or CrsPick==4:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip=-1)
							elif CrsPick==5 or CrsPick==6:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 1)
							else:
								FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],
								                     SubC=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 0)
							FixVCrsNet.RandomRest()
							FixVCrsNet.SetName(Name=Name(Prefix='FixVCrsNet', Index=SampleIndex))
							FixVCrsNet.Write(Folder+'FixVCrsNet/')
							SampleIndex+=1
							# ---------------------------------------------------------------------
	print(SampleIndex)

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
# CreatePoreNetworkSamples(Nx= 3, Ny= 3, Folder='/home/xu/work/PoreNetwork1010Samples')
CreatePoreNetworkSamples(Nx=10, Ny=10, Folder='/home/xu/work/PoreNetwork1010Samples/')
CreatePoreNetworkSamples(Nx=20, Ny=20, Folder='/home/xu/work/PoreNetwork2020Samples/')
CreatePoreNetworkSamples(Nx=40, Ny=40, Folder='/home/xu/work/PoreNetwork4040Samples/')