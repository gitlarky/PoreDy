#==================================================================================================
# Title      : SmartNet
# Description: Pore-Network Database Generation
# Author     : Zhenyu Xu; westlark@outlook.com
# Start Time : 2019.09.01
# License    : Apache 2.0
#==================================================================================================

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

# List Operation ----------------------------------------------------------------------------------
def SetUnion(Base=[], Sub=[]):
	final=[]
	for item in Base:
		final.append(item)
	for item in Sub:
		if not item in Base:
			final.append(item)
	return final

# List Operation ----------------------------------------------------------------------------------
def SetJoint(Base=[], Sub=[]):
	final=[]
	for item in Base:
		if item in Sub:
			final.append(item)
	return final

# List Operation ----------------------------------------------------------------------------------
def SetDifference(Base=[], Exc=[]):
	final=[]
	for item in Base:
		if item in Exc:
			continue
		else:
			final.append(item)
	return final



#============================ Class ThroatPool ====================================================
class ThroatPool(object):
	# Initialize ----------------------------------------------------------------------------------
	def __init__(self, Name='ThroatPool', ThroatType=[], ThroatCount=[]):
		self.Name        =Name
		self.TypeCount   =min(len(ThroatType), len(ThroatCount))
		self.Pool        ={}
		self.InitialPool ={}
		self.InitialTotal=0
		for i in range(self.TypeCount):
			self.Pool[ThroatType[i]]=ThroatCount[i]
			self.InitialPool[ThroatType[i]]=ThroatCount[i]
			self.InitialTotal+=ThroatCount[i]

	# Pool Total Throat Count ---------------------------------------------------------------------
	def ExistTotal(self):
		et=0
		for count in self.Pool.values():
			et+=count
		return et

	# Throat used count ---------------------------------------------------------------------------
	def UsedTotal(self):
		return self.InitialTotal-self.ExistTotal()

	# Find out the index of existing throat -------------------------------------------------------
	def ExistType(self):
		et=[]
		for throat, count in self.Pool.items():
			if count>0:
				et.append(throat)
		return et

	# Find out the existing throat pool -------------------------------------------------------
	def ExistPool(self):
		return {k:v for k, v in self.Pool.items() if v>0}

	# Report the pool -----------------------------------------------------------------------------
	def Report(self):
		print(self.Name, ':', self.InitialTotal, self.ExistTotal(), self.UsedTotal()\
			           , '\t', self.InitialPool\
			           , '\t', self.Pool, self.ExistType())
		return True

	# Throat in Pool or not -----------------------------------------------------------------------
	def InPool(self, T):
		T=abs(T)
		if self.Pool.get(T)!=None:
			return True
		else:
			return False

	# Fetch one -----------------------------------------------------------------------------------
	def FetchOne(self, Index=0, Sub=[], Exc=[], Method='Random'):
		Base=self.ExistType()
		if not Base:
			print('Warning: -------------------------------------------------Throat Pool Used Up!')
			return 0
		else:
			if Method=='Random':
				Selection=Pick(List=Base)
			elif Method=='Rotate':
				if not Sub and not Exc:
					Selection=Pick(List=Base, Index=Index, Method='Rotate')
				elif Sub and not Exc:
					Selection=Pick(List=Sub, Index=Index, Method='Rotate')
					if not Selection in Base:
						Selection=Pick(List=Base)
				elif not Sub and Exc:
					Selection=Pick(List=Base, Index=Index, Method='Rotate')
					if Selection in Exc:
						diff=SetDifference(Base, Exc)
						if diff:
							Selection=Pick(List=diff)
				elif Sub and Exc:
					Selection=Pick(List=Sub, Index=Index, Method='Rotate')
					if not Selection in Base:
						diff=SetDifference(Base, Exc)
						if diff:
							Selection=Pick(List=diff)
						else:
							Selection=Pick(List=Base)
			self.Pool[Selection]-=1
			return Selection

	# Return One ----------------------------------------------------------------------------------
	def ReturnOne(self, T):
		T=abs(T)
		if self.InPool(T):
			self.Pool[T]+=1
			return True
		else:
			print('Warning: cannot return, no throat type matched!')
			return False

#============================ Class PoreNetwork ===================================================
class PoreNetwork(object):
	# Initialize ----------------------------------------------------------------------------------
	def __init__(self, Name='PoreNetwork', Nx=20, Ny=20, Open=['N'], Cross=False, FixV=False, 
		         StrTType=[], CrsTType=[]):
		self.Name     =Name
		self.Nx       =Nx
		self.Ny       =Ny
		self.Open     =Open
		self.Cross    =Cross
		self.FixV     =FixV

		self.VT       =np.zeros((self.Ny+1, self.Nx+1), dtype=int)
		self.HT       =np.zeros((self.Ny+1, self.Nx+1), dtype=int)
		self.CT       =np.zeros((self.Ny+1, self.Nx+1), dtype=int)
		self.VTRange  =[[0, self.Nx+1], [1, self.Ny+1]]
		self.HTRange  =[[1, self.Nx+1], [0, self.Ny+1]]
		self.CTRange  =[[1, self.Nx+1], [1, self.Ny+1]] if self.Cross else [[0, 0], [0, 0]]
		self.StrTCount=self.Ny*(self.Nx+1)+self.Nx*(self.Ny+1)
		self.CrsTCount=self.Nx*self.Ny if self.Cross else 0
		self.TotTCount=self.StrTCount+self.CrsTCount
		if 'W' in self.Open:
			self.VTRange[0][0] =1
			self.StrTCount    -=self.Ny
		if 'E' in self.Open:
			self.VTRange[0][1] =self.Nx-1
			self.StrTCount    -=self.Ny
		if 'S' in self.Open:
			self.HTRange[1][0] =1
			self.StrTCount    -=self.Nx
		if 'N' in self.Open:
			self.HTRange[1][1] =self.Ny-1
			self.StrTCount    -=self.Nx

		STN=[self.StrTCount for i in range(len(StrTType))]
		if self.FixV:
			AvgStrTN              =self.StrTCount//len(StrTType)
			RemStrTN              =self.StrTCount%len(StrTType)
			STN                   =[AvgStrTN for i in range(len(StrTType))]
			STN[len(StrTType)//2]+=RemStrTN
		self.StrTPool=ThroatPool(Name=self.Name+'.StrTPool', ThroatType=StrTType, ThroatCount=STN)
		if self.Cross:
			CTN=[self.CrsTCount for i in range(len(CrsTType))]
			if self.FixV:
				AvgCrsTN              =self.CrsTCount//len(CrsTType)
				RemCrsTN              =self.CrsTCount%len(CrsTType)
				CTN                   =[AvgCrsTN for i in range(len(CrsTType))]
				CTN[len(CrsTType)//2]+=RemCrsTN
			self.CrsTPool=ThroatPool(Name=self.Name+'.CrsTPool', ThroatType=CrsTType, ThroatCount=CTN)

		if self.StrTCount>self.StrTPool.ExistTotal():
			print('Warning: Straight Throat Pool Created but Not Enough!')
		if self.Cross and self.CrsTCount>self.CrsTPool.ExistTotal():
			print('Warning: Cross Throat Pool Created but Not Enough!')

	# set a new name ------------------------------------------------------------------------------
	def SetName(self, Name):
		self.Name=Name
		return True

	# Output the file binary and other related files ----------------------------------------------
	def Output(self, Folder='', Disp=False, Dump=True, Write=False, Pixel=False, ColorMap='Gray'):
		if Dump: pickle.dump(self, open(Folder+self.Name+'.net', 'wb'))
		if Disp: print('Successfully saved pore-network binary file.')

		if Write or Pixel:
			Mx=2*self.Nx+1
			My=2*self.Ny+1
			Matrix=np.zeros((My, Mx), dtype=int)
			for J in range(self.VTRange[1][0], self.VTRange[1][1], 1):
				for I in range(self.VTRange[0][0], self.VTRange[0][1], 1):
					Matrix[J*2-1][I*2  ]=self.VT[J][I]
			for J in range(self.HTRange[1][0], self.HTRange[1][1], 1):
				for I in range(self.HTRange[0][0], self.HTRange[0][1], 1):
					Matrix[J*2  ][I*2-1]=self.HT[J][I]
			for J in range(self.CTRange[1][0], self.CTRange[1][1], 1):
				for I in range(self.CTRange[0][0], self.CTRange[0][1], 1):
					Matrix[J*2-1][I*2-1]=self.CT[J][I]

			if Write:
				with open(Folder+self.Name+'.at', 'w') as wat:
					for j in range(self.My-1, -1, -1):
						for i in range(self.Mx):
							if Matrix[j][i]==0:
								wat.write('% 9.6e\t% 9d\t' % (0, 0))
							else:
								throattype=Matrix[j][i]
								flip=throattype/abs(throattype)
								row =throattype%10
								col =throattype//10
								wat.write('% 9.6e\t% 9d\t' % (ThroatChoice[row][col].D*flip, ThroatChoice[row][col].S))
						wat.write('\n')

			if Pixel:
				if ColorMap=='Gray':
					D=np.zeros((self.My, self.Mx), dtype=float)
					for j in range(self.My-1, -1, -1):
						for i in range(self.Mx):
							throattype=Matrix[j][i]
							flip      =throattype/abs(throattype)
							row       =throattype%10
							col       =throattype//10
							D[j][i]   =ThroatChoice[row][col].D*flip
					plt.imshow(D, cmap="gray", vmin=0, vmax=225)
					if Disp: plt.show()
					if Save: plt.savefig(Folder+self.Name+'.png')
				elif ColorMap=='RGB':
					pass

		return True

	# Get Throat info at a certain position -------------------------------------------------------
	def GetT(self, TP, I, J):
		if TP=='VT':
			if (I>=self.VTRange[0][0] and I<self.VTRange[0][1] and J>=self.VTRange[1][0] and J<self.VTRange[1][1]):
				return ['VT', I, J, self.VT[J][I]]
			else:
				return ['VT', I, J]
		elif TP=='HT':
			if (I>=self.HTRange[0][0] and I<self.HTRange[0][1] and J>=self.HTRange[1][0] and J<self.HTRange[1][1]):
				return ['HT', I, J, self.HT[J][I]]
			else:
				return ['HT', I, J]
		elif TP=='CT':
			if (I>=self.CTRange[0][0] and I<self.CTRange[0][1] and J>=self.CTRange[1][0] and J<self.CTRange[1][1]) and self.Cross:
				return ['CT', I, J, self.CT[J][I]]
			else:
				return ['CT', I, J]
		else:
			return []

	# Assign Throat at a certain position ---------------------------------------------------------
	def AssignT(self, TP, I, J, T):
		if   TP=='VT' and I>=self.VTRange[0][0] and I<self.VTRange[0][1] and J>=self.VTRange[1][0] and J<self.VTRange[1][1]:
			self.VT[J][I]=T
		elif TP=='HT' and I>=self.HTRange[0][0] and I<self.HTRange[0][1] and J>=self.HTRange[1][0] and J<self.HTRange[1][1]:
			self.HT[J][I]=T
		elif TP=='CT' and I>=self.CTRange[0][0] and I<self.CTRange[0][1] and J>=self.CTRange[1][0] and J<self.CTRange[1][1] and self.Cross:
			self.CT[J][I]=T
		else:
			return False
		return True

	# Judge if it is assigned or not --------------------------------------------------------------
	def Assigned(self, TP, I, J):
		thro=self.GetT(TP, I, J)
		if len(thro)==4:
			if TP=='VT' or TP=='HT':
				if self.StrTPool.InPool(thro[3]):
					return True
				else:
					return False
			elif TP=='CT':
				if self.Cross:
					if self.CrsTPool.InPool(thro[3]):
						return True
					else:
						return False
				else:
					return True
		else:
			print('Warning: There should be no Throat, Donot Operate Me!')
			return True
	
	# Return Throat at a certain position back to appropriate pool --------------------------------
	def ReturnT(self, TP, I, J):
		T=self.GetT(TP, I, J)
		if  (TP=='VT' or TP=='HT') and len(T)==4:
			if T[3]==0: print('Warning: where is this throat coming from?: ', TP, I, J, T)
			succeed=self.StrTPool.ReturnOne(T=T[3])
			if not succeed: print('Warning: Returning Operation Failed!')
			self.AssignT(TP=TP, I=I, J=J, T=0)
		elif TP=='CT' and len(T)==4:
			if T[3]==0: print('Warning: where is this throat coming from?: ', TP, I, J, T)
			succeed=self.CrsTPool.ReturnOne(T=T[3])
			if not succeed: print('Warning: Returning Operation Failed!')
			self.AssignT(TP=TP, I=I, J=J, T=0)
		else:
			print('Warning: No Throat to Return!')
			return False
		return True

	# Restore all throats back to appropriate pool ------------------------------------------------
	def Restore(self, TPs=['VT', 'HT', 'CT']):
		vt=0
		ht=0
		ct=0
		if 'VT' in TPs:
			for I in range(self.VTRange[0][0], self.VTRange[0][1], 1):
				for J in range(self.VTRange[1][0], self.VTRange[1][1], 1):
					if self.ReturnT('VT', I, J): vt+=1
		if 'HT' in TPs:
			for I in range(self.HTRange[0][0], self.HTRange[0][1], 1):
				for J in range(self.HTRange[1][0], self.HTRange[1][1], 1):
					if self.ReturnT('HT', I, J): ht+=1
		if 'CT' in TPs and self.Cross:
			for I in range(self.CTRange[0][0], self.CTRange[0][1], 1):
				for J in range(self.CTRange[1][0], self.CTRange[1][1], 1):
					if self.ReturnT('CT', I, J): ct+=1
		print('Restore: ', vt+ht+ct, '\tVT: ', vt, '\tHT: ', ht, '\tCT: ', ct)
		return True

	# Assign Throat distribution in a box of certain position -------------------------------------
	def AssignBox(self, TP='',
		                Start=[0, 0], End=[0, 0], Band=[0, 0],
		                Sub=[], Grad=[0, 0, 0], Repeat=[1, 1, 1], Flip=0):
		if   TP=='VT':
			Range=self.VTRange
			PoolChoice=self.StrTPool
		elif TP=='HT':
			Range=self.HTRange
			PoolChoice=self.StrTPool
		elif TP=='CT' and self.Cross:
			Range=self.CTRange
			PoolChoice=self.CrsTPool
		if Start[0]==0 and End[0]==0:
			Start[0]=Range[0][0]
			End  [0]=Range[0][1]
		if Start[1]==0 and End[1]==0:
			Start[1]=Range[1][0]
			End  [1]=Range[1][1]
		if Sub==[]:
			Sub=[key for key in PoolChoice.Pool.keys()]
			Exc=[]
		else:
			Exc=Sub

		Count=0
		Replace=0
		for I in range(Start[0], End[0], 1):
			for J in range(Start[1], End[1], 1):
				if Band==[0, 0] or (\
				   (Band[0]> 0 and Band[1]> 0) and not(I>=Start[0]+Band[0] and I<=End[0]-Band[0]-1 and \
				                                       J>=Start[1]+Band[1] and J<=End[1]-Band[1]-1)) or (\
				   (Band[0]> 0 and Band[1]==0) and not(I>=Start[0]+Band[0] and I<=End[0]-Band[0]-1)) or (\
				   (Band[0]==0 and Band[1]> 0) and not(J>=Start[1]+Band[1] and J<=End[1]-Band[1]-1)):
					if self.Assigned(TP, I, J):
						self.ReturnT(TP, I, J)
						Replace+=1

					for ring in range(min(Band[0], Band[1])+1):
						if not (I>=Start[0]+ring+1 and I<=End[0]-ring-2 and 
							J>=Start[1]+ring+1 and J<=End[1]-ring-2) and (
							I>=Start[0]+ring   and I<=End[0]-ring-1 and 
							J>=Start[1]+ring   and J<=End[1]-ring-1):
							RingIndex=ring
					Selection=Pick(List=Sub, Index=(I-Start[0])//Repeat[0]*Grad[0]+\
						                            (J-Start[1])//Repeat[1]*Grad[1]+\
						                             RingIndex  //Repeat[2]*Grad[2], Method='Rotate')
					if TP=='VT' or TP=='HT':
						T=self.StrTPool.FetchOne(Index=Selection, Sub=Sub, Exc=Exc, Method='Rotate')
					elif TP=='CT':
						T=self.CrsTPool.FetchOne(Index=Selection, Sub=Sub, Exc=Exc, Method='Rotate')
						if not self.CrsTPool.InPool(T): print('Warning: Fuck it is impossible!')
						if Flip==0:
							T= T*Pick(List=[-1, 1])
						elif Flip==2 and I%2==0:
							T=-T
						elif Flip==3 and J%2==0:
							T=-T
						elif Flip==-1 or Flip==1: # Flip= -1 or 1
							T= T*Flip
						if not self.CrsTPool.InPool(T): print('Warning: Fuck! What happened here!', t, T)
					succeed=self.AssignT(TP, I, J, T)
					if not succeed: print('Warning: AssignT when AssignBox not succeed!')
					Count+=1
		print('AssignBox: ', TP, Start, End, '\tCount=', Count, '\tReplace=', Replace)
		return Count

	# Slice a Region from original Matrix ---------------------------------------------------------
	def RandomRest(self):
		vt=0
		ht=0
		ct=0
		for I in range(self.VTRange[0][0], self.VTRange[0][1], 1):
			for J in range(self.VTRange[1][0], self.VTRange[1][1], 1):
				if not self.Assigned('VT', I, J):
					succeed=self.AssignT('VT', I, J, self.StrTPool.FetchOne())
					if succeed:
						vt+=1
					else:
						print('Warning: AssignT when RandomRest not succeed!')
		for I in range(self.HTRange[0][0], self.HTRange[0][1], 1):
			for J in range(self.HTRange[1][0], self.HTRange[1][1], 1):
				if not self.Assigned('HT', I, J):
					succeed=self.AssignT('HT', I, J, self.StrTPool.FetchOne())
					if succeed:
						ht+=1
					else:
						print('Warning: AssignT when RandomRest not succeed!')
		if self.Cross:
			for I in range(self.CTRange[0][0], self.CTRange[0][1], 1):
				for J in range(self.CTRange[1][0], self.CTRange[1][1], 1):
					if not self.Assigned('CT', I, J):
						succeed=self.AssignT('CT', I, J, self.CrsTPool.FetchOne()*Pick(List=[-1, 1]))
						if succeed:
							ct+=1
						else:
							print('Warning: AssignT when RandomRest not succeed!')
		print('RandomRest: ', vt+ht+ct, '\tVT: ', vt, '\tHT: ', ht, '\tCT: ', ct)
		return vt+ht+ct

	# Check if it is correct ----------------------------------------------------------------------
	def Check(self):
		vt=0
		ht=0
		ct=0
		for TP in ['HT', 'VT', 'CT']:
			if TP=='VT':
				for I in range(self.VTRange[0][0], self.VTRange[0][1], 1):
					for J in range(self.VTRange[1][0], self.VTRange[1][1], 1):
						throat=self.GetT(TP=TP, I=I, J=J)
						if not self.StrTPool.InPool(throat[3]):
							print('Error: throat not correctly assigned!---------------------------', TP, I, J, throat)
							return False
			if TP=='HT':
				for I in range(self.HTRange[0][0], self.HTRange[0][1], 1):
					for J in range(self.HTRange[1][0], self.HTRange[1][1], 1):
						throat=self.GetT(TP=TP, I=I, J=J)
						if not self.StrTPool.InPool(throat[3]):
							print('Error: throat not correctly assigned!---------------------------', TP, I, J, throat)
							return False
			if TP=='CT' and self.Cross:
				for I in range(self.CTRange[0][0], self.CTRange[0][1], 1):
					for J in range(self.CTRange[1][0], self.CTRange[1][1], 1):
						throat=self.GetT(TP=TP, I=I, J=J)
						if not self.StrTPool.InPool(throat[3]):
							print('Error: throat not correctly assigned!---------------------------', TP, I, J, throat)
							return False
		return True

	# Report information --------------------------------------------------------------------------
	def Report(self):
		print(self.Name, ':', '\tNx=', self.Nx, 'Ny=', self.Ny, '\tTotTCount=', self.TotTCount, '\tFixV=', self.FixV)
		print('\tOpen: ', self.Open , '\tStrTCount=', self.StrTCount, '\tVTRange=', self.VTRange, 'HTRange=', self.HTRange)
		print('\tStrTType= ', self.StrTType, '\tStrTPool: ', self.StrTPool.Pool)
		print('\tCross=', self.Cross, '\tCrsTCount=', self.CrsTCount, '\tCTRange=', self.CTRange, )
		print('\tCrsTType= ', self.CrsTType, '\tCrsTPool: ', self.CrsTPool.Pool)

		return True		

#============================ Create Pore-Network Samples =========================================
# Create samples for One Assignment Region for each of the 4 kinds of networks ================
def CreateOne(Nx=20, Ny=20, Folder='', SIN=0, PNW=['Rand', 'FixV'], AVTP=0, AHTP=0, \
	          IS=0, JS=0, IE=0, JE=0, IB=0, JB=0, \
	          SC=[], IG=0, JG=0, RG=0, IR=1, JR=1, RR=1):
	if 'Rand' in PNW:
		# -------------------------------------------------------------------------------------
		RandStrNet =PoreNetwork(Name=Name(Prefix='RandStrNet', Index=SIN),\
		                        Nx=Nx, Ny=Ny, Open=['N'], Cross=False, FixV=False,\
		                        StrTType=TT)
		if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<AVTP:
			RandStrNet.AssignBox(TP='VT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],\
				                 Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
		if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<AHTP:
			RandStrNet.AssignBox(TP='HT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],\
				                 Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
		RandStrNet.RandomRest()
		RandStrNet.Output(Folder+'RandStrNet/')
		RandStrNet.Check()
		print('RandStrNet', SIN, 'Written!:\t', [IS, IE], [JS, JE], [IB, JB], SC, [IG, JG, RG], [IR, JR, RR])
		# -------------------------------------------------------------------------------------
		RandCrsNet =PoreNetwork(Name=Name(Prefix='RandCrsNet', Index=SIN), \
		                        Nx=Nx, Ny=Ny, Open=['N'], Cross=True , FixV=False, \
		                        StrTType=TT, CrsTType=TT)
		if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<AVTP:
			RandCrsNet.AssignBox(TP='VT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],\
				                 Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
		if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<AHTP:
			RandCrsNet.AssignBox(TP='HT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],\
				                 Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
		CrsPick=Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
		if CrsPick<2:
			RandCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 2)
		elif CrsPick>7:
			RandCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 3)
		elif CrsPick==3:
			RandCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip=-1)
		elif CrsPick==6:
			RandCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 1)
		elif CrsPick==4:
			RandCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],\
			                     Sub=[], Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip=-1)
		elif CrsPick==5:
			RandCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],\
			                     Sub=[], Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 1)
		else:
			RandCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 0)
		RandCrsNet.RandomRest()
		RandCrsNet.Output(Folder+'RandCrsNet/')
		RandCrsNet.Check()
		print('RandCrsNet', SIN, 'Written!:\t', [IS, IE], [JS, JE], [IB, JB], SC, [IG, JG, RG], [IR, JR, RR])
		# -------------------------------------------------------------------------------------
	if 'FixV' in PNW:
		# -------------------------------------------------------------------------------------
		FixVStrNet =PoreNetwork(Name=Name(Prefix='FixVStrNet', Index=SIN),\
		                        Nx=Nx, Ny=Ny, Open=['N'], Cross=False, FixV=True,\
		                        StrTType=TT)
		if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<AVTP:
			FixVStrNet.AssignBox(TP='VT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],\
				                 Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
		if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<AHTP:
			FixVStrNet.AssignBox(TP='HT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],\
				                 Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
		FixVStrNet.RandomRest()
		FixVStrNet.Output(Folder+'FixVStrNet/')
		FixVStrNet.Check()
		print('FixVStrNet', SIN, 'Written!:\t', [IS, IE], [JS, JE], [IB, JB], SC, [IG, JG, RG], [IR, JR, RR])
		# -------------------------------------------------------------------------------------
		FixVCrsNet =PoreNetwork(Name=Name(Prefix='FixVCrsNet', Index=SIN), \
		                        Nx=Nx, Ny=Ny, Open=['N'], Cross=True , FixV=True, \
		                        StrTType=TT, CrsTType=TT)
		if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<AVTP:
			FixVCrsNet.AssignBox(TP='VT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],\
				                 Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
		if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<AHTP:
			FixVCrsNet.AssignBox(TP='HT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],\
				                 Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
		CrsPick=Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
		if CrsPick<2:
			FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 2)
		elif CrsPick>7:
			FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 3)
		elif CrsPick==3:
			FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip=-1)
		elif CrsPick==6:
			FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 1)
		elif CrsPick==4:
			FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],\
			                     Sub=[], Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip=-1)
		elif CrsPick==5:
			FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],\
			                     Sub=[], Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 1)
		else:
			FixVCrsNet.AssignBox(TP='CT', Start=[IS, JS], End=[IE, JE], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 0)
		FixVCrsNet.RandomRest()
		FixVCrsNet.Output(Folder+'FixVCrsNet/')
		FixVCrsNet.Check()
		print('FixVCrsNet', SIN, 'Written!:\t', [IS, IE], [JS, JE], [IB, JB], SC, [IG, JG, RG], [IR, JR, RR])
		# -------------------------------------------------------------------------------------
	return True
# =============================================================================================
# Create samples for Two Assignment Region for each of the 4 kinds of networks ================
def CreateTwo(Nx=20, Ny=20, Folder='', SIN=0, PNW=['Rand', 'FixV'], AVTP=0, AHTP=0, \
	          IS=0, JS=0, IE=0, JE=0, IB=0, JB=0, \
	          SC=[], IG=0, JG=0, RG=0, IR=1, JR=1, RR=1):
	if 'Rand' in PNW:
		# -------------------------------------------------------------------------------------
		RandStrNet =PoreNetwork(Name=Name(Prefix='RandStrNet', Index=SIN),\
		                        Nx=Nx, Ny=Ny, Open=['N'], Cross=False, FixV=False,\
		                        StrTType=TT)
		if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<AVTP:
			RandStrNet.AssignBox(TP='VT', Start=[IS, 0], End=[IE, 0], Band=[IB, JB],\
				                 Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
			RandStrNet.AssignBox(TP='VT', Start=[0, JS], End=[0, JE], Band=[IB, JB],\
				                 Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
		if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<AHTP:
			RandStrNet.AssignBox(TP='HT', Start=[0, JS], End=[0, JE], Band=[IB, JB],\
				                 Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
			RandStrNet.AssignBox(TP='HT', Start=[IS, 0], End=[IE, 0], Band=[IB, JB],\
				                 Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
		RandStrNet.RandomRest()
		RandStrNet.Output(Folder+'RandStrNet/')
		RandStrNet.Check()
		print('RandStrNet', SIN, 'Written!:\t', [IS, IE], [JS, JE], [IB, JB], SC, [IG, JG, RG], [IR, JR, RR])
		# -------------------------------------------------------------------------------------
		RandCrsNet =PoreNetwork(Name=Name(Prefix='RandCrsNet', Index=SIN), \
		                        Nx=Nx, Ny=Ny, Open=['N'], Cross=True , FixV=False, \
		                        StrTType=TT, CrsTType=TT)
		if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<AVTP:
			RandCrsNet.AssignBox(TP='VT', Start=[0, JS], End=[0, JE], Band=[IB, JB],\
				                 Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
			RandCrsNet.AssignBox(TP='VT', Start=[IS, 0], End=[IE, 0], Band=[IB, JB],\
				                 Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
		if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<AHTP:
			RandCrsNet.AssignBox(TP='HT', Start=[IS, 0], End=[IE, 0], Band=[IB, JB],\
				                 Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
			RandCrsNet.AssignBox(TP='HT', Start=[0, JS], End=[0, JE], Band=[IB, JB],\
				                 Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
		CrsPick=Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
		if CrsPick<2:
			RandCrsNet.AssignBox(TP='CT', Start=[0, JS], End=[0, JE], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 2)
			RandCrsNet.AssignBox(TP='CT', Start=[IS, 0], End=[IE, 0], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 2)
		elif CrsPick>7:
			RandCrsNet.AssignBox(TP='CT', Start=[IS, 0], End=[IE, 0], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 3)
			RandCrsNet.AssignBox(TP='CT', Start=[0, JS], End=[0, JE], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 3)
		elif CrsPick==3:
			RandCrsNet.AssignBox(TP='CT', Start=[0, JS], End=[0, JE], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip=-1)
			RandCrsNet.AssignBox(TP='CT', Start=[IS, 0], End=[IE, 0], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip=-1)
		elif CrsPick==6:
			RandCrsNet.AssignBox(TP='CT', Start=[IS, 0], End=[IE, 0], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 1)
			RandCrsNet.AssignBox(TP='CT', Start=[0, JS], End=[0, JE], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 1)
		elif CrsPick==4:
			RandCrsNet.AssignBox(TP='CT', Start=[0, JS], End=[0, JE], Band=[IB, JB],\
			                     Sub=[], Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip=-1)
			RandCrsNet.AssignBox(TP='CT', Start=[IS, 0], End=[IE, 0], Band=[IB, JB],\
			                     Sub=[], Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip=-1)
		elif CrsPick==5:
			RandCrsNet.AssignBox(TP='CT', Start=[IS, 0], End=[IE, 0], Band=[IB, JB],\
			                     Sub=[], Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 1)
			RandCrsNet.AssignBox(TP='CT', Start=[0, JS], End=[0, JE], Band=[IB, JB],\
			                     Sub=[], Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 1)
		else:
			RandCrsNet.AssignBox(TP='CT', Start=[0, JS], End=[0, JE], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 0)
			RandCrsNet.AssignBox(TP='CT', Start=[IS, 0], End=[IE, 0], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 0)
		RandCrsNet.RandomRest()
		RandCrsNet.Output(Folder+'RandCrsNet/')
		RandCrsNet.Check()
		print('RandCrsNet', SIN, 'Written!:\t', [IS, IE], [JS, JE], [IB, JB], SC, [IG, JG, RG], [IR, JR, RR])
		# -------------------------------------------------------------------------------------
	if 'FixV' in PNW:
		# -------------------------------------------------------------------------------------
		FixVStrNet =PoreNetwork(Name=Name(Prefix='FixVStrNet', Index=SIN),\
		                        Nx=Nx, Ny=Ny, Open=['N'], Cross=False, FixV=True,\
		                        StrTType=TT)
		if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<AVTP:
			FixVStrNet.AssignBox(TP='VT', Start=[0, JS], End=[0, JE], Band=[IB, JB],\
				                 Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
			FixVStrNet.AssignBox(TP='VT', Start=[IS, 0], End=[IE, 0], Band=[IB, JB],\
				                 Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
		if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<AHTP:
			FixVStrNet.AssignBox(TP='HT', Start=[IS, 0], End=[IE, 0], Band=[IB, JB],\
				                 Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
			FixVStrNet.AssignBox(TP='HT', Start=[0, JS], End=[0, JE], Band=[IB, JB],\
				                 Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
		FixVStrNet.RandomRest()
		FixVStrNet.Output(Folder+'FixVStrNet/')
		FixVStrNet.Check()
		print('FixVStrNet', SIN, 'Written!:\t', [IS, IE], [JS, JE], [IB, JB], SC, [IG, JG, RG], [IR, JR, RR])
		# -------------------------------------------------------------------------------------
		FixVCrsNet =PoreNetwork(Name=Name(Prefix='FixVCrsNet', Index=SIN), \
		                        Nx=Nx, Ny=Ny, Open=['N'], Cross=True , FixV=True, \
		                        StrTType=TT, CrsTType=TT)
		if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<AVTP:
			FixVCrsNet.AssignBox(TP='VT', Start=[IS, 0], End=[IE, 0], Band=[IB, JB],\
				                 Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
			FixVCrsNet.AssignBox(TP='VT', Start=[0, JS], End=[0, JE], Band=[IB, JB],\
				                 Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
		if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<AHTP:
			FixVCrsNet.AssignBox(TP='HT', Start=[0, JS], End=[0, JE], Band=[IB, JB],\
				                 Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
			FixVCrsNet.AssignBox(TP='HT', Start=[IS, 0], End=[IE, 0], Band=[IB, JB],\
				                 Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR])
		CrsPick=Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
		if CrsPick<2:
			FixVCrsNet.AssignBox(TP='CT', Start=[IS, 0], End=[IE, 0], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 2)
			FixVCrsNet.AssignBox(TP='CT', Start=[0, JS], End=[0, JE], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 2)
		elif CrsPick>7:
			FixVCrsNet.AssignBox(TP='CT', Start=[0, JS], End=[0, JE], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 3)
			FixVCrsNet.AssignBox(TP='CT', Start=[IS, 0], End=[IE, 0], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 3)
		elif CrsPick==3:
			FixVCrsNet.AssignBox(TP='CT', Start=[IS, 0], End=[IE, 0], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip=-1)
			FixVCrsNet.AssignBox(TP='CT', Start=[0, JS], End=[0, JE], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip=-1)
		elif CrsPick==6:
			FixVCrsNet.AssignBox(TP='CT', Start=[0, JS], End=[0, JE], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 1)
			FixVCrsNet.AssignBox(TP='CT', Start=[IS, 0], End=[IE, 0], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 1)
		elif CrsPick==4:
			FixVCrsNet.AssignBox(TP='CT', Start=[IS, 0], End=[IE, 0], Band=[IB, JB],\
			                     Sub=[], Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip=-1)
			FixVCrsNet.AssignBox(TP='CT', Start=[0, JS], End=[0, JE], Band=[IB, JB],\
			                     Sub=[], Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip=-1)
		elif CrsPick==5:
			FixVCrsNet.AssignBox(TP='CT', Start=[0, JS], End=[0, JE], Band=[IB, JB],\
			                     Sub=[], Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 1)
			FixVCrsNet.AssignBox(TP='CT', Start=[IS, 0], End=[IE, 0], Band=[IB, JB],\
			                     Sub=[], Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 1)
		else:
			FixVCrsNet.AssignBox(TP='CT', Start=[IS, 0], End=[IE, 0], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 0)
			FixVCrsNet.AssignBox(TP='CT', Start=[0, JS], End=[0, JE], Band=[IB, JB],\
			                     Sub=SC, Grad=[IG, JG, RG], Repeat=[IR, JR, RR], Flip= 0)
		FixVCrsNet.RandomRest()
		FixVCrsNet.Output(Folder+'FixVCrsNet/')
		RandCrsNet.Check()
		print('FixVCrsNet', SIN, 'Written!:\t', [IS, IE], [JS, JE], [IB, JB], SC, [IG, JG, RG], [IR, JR, RR])
		# -------------------------------------------------------------------------------------
	return True
# =============================================================================================
def CreatePoreNetworkSamples(Nx=20, Ny=20, Folder=''):
	print('Prepare the folders ------------------------------------------------------------------')
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

	print('Generate Pore-Networks ---------------------------------------------------------------')
	# SampleIndex=0
	# print('Generating 1st Group of Pore-Networks ------------------------------------------------')
	# for SC in [[0, 1], [2, 3], [6, 5], [8, 7]]:
	# 	for IS in range(1, Nx-7, 4):
	# 		IE=IS+8
	# 		for IG in [0, 1]:
	# 			for IR in [1, 3]:
	# 				for JR in [2, 4]:
	# 					CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IS=IS, IE=IE, SC=SC, IG=IG, IR=IR, JR=JR, AVTP=8)
	# 					SampleIndex+=1
	# 					CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IS=IS, IE=IE, SC=SC, IG=IG, IR=IR, JR=JR, AHTP=8)
	# 					SampleIndex+=1
	# 					CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IS=IS, IE=IE, SC=SC, IG=IG, IR=IR, JR=JR, AVTP=9, AHTP=9)
	# 					SampleIndex+=1
	# 	for JS in range(1, Ny-7, 4):
	# 		JE=JS+8
	# 		for JG in [0, 1]:
	# 			for IR in [2, 3]:
	# 				for JR in [1, 3]:
	# 					CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, JS=JS, JE=JE, SC=SC, JG=JG, IR=IR, JR=JR, AVTP=8)
	# 					SampleIndex+=1
	# 					CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, JS=JS, JE=JE, SC=SC, JG=JG, IR=IR, JR=JR, AHTP=8)
	# 					SampleIndex+=1
	# 					CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, JS=JS, JE=JE, SC=SC, JG=JG, IR=IR, JR=JR, AVTP=9, AHTP=9)
	# 					SampleIndex+=1
	# 	for IS in range(1, Nx-7, 4):
	# 		IE=IS+8
	# 		for JS in range(1, Ny-7, 4):
	# 			JE=JS+8
	# 			IG=1
	# 			JG=1
	# 			for IR in [1, 2, 3, 4]:
	# 				for JR in [1, 2, 3, 4]:
	# 					CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AVTP=8)
	# 					SampleIndex+=1
	# 					CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AHTP=8)
	# 					SampleIndex+=1
	# 					CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AVTP=9, AHTP=9)
	# 					SampleIndex+=1
	# 	for IB in [1, 2, 3, 4]:
	# 		for IG in [0, 1]:
	# 			for IR in [2, 4]:
	# 				for JR in [1, 3]:
	# 					CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IB=IB, SC=SC, IG=IG, IR=IR, JR=JR, AVTP=8)
	# 					SampleIndex+=1
	# 					CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IB=IB, SC=SC, IG=IG, IR=IR, JR=JR, AHTP=8)
	# 					SampleIndex+=1
	# 					CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IB=IB, SC=SC, IG=IG, IR=IR, JR=JR, AVTP=9, AHTP=9)
	# 					SampleIndex+=1
	# 	for JB in [1, 2, 3, 4]:
	# 		for JG in [0, 1]:
	# 			for IR in [1, 3]:
	# 				for JR in [2, 4]:
	# 					CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, JB=JB, SC=SC, JG=JG, IR=IR, JR=JR, AVTP=8)
	# 					SampleIndex+=1
	# 					CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, JB=JB, SC=SC, JG=JG, IR=IR, JR=JR, AHTP=8)
	# 					SampleIndex+=1
	# 					CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, JB=JB, SC=SC, JG=JG, IR=IR, JR=JR, AVTP=9, AHTP=9)
	# 					SampleIndex+=1
	# print('Group 1 of Pore-Networks Generated: ', SampleIndex, '---------------------------------')
	SampleIndex=3072
	print('Generating 2nd Group of Pore-Networks ------------------------------------------------')
	for SC in [[0, 1, 2], [4, 3, 5], [8, 7, 6]]:
		for IS in range(1, Nx-11, 4):
			IE=IS+12
			for JG in [0, 1, 2]:
				for IR in [2, 3, 4]:
					for JR in [2, 3, 4]:
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IS=IS, IE=IE, SC=SC, JG=JG, IR=IR, JR=JR, AVTP=8)
						SampleIndex+=1
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IS=IS, IE=IE, SC=SC, JG=JG, IR=IR, JR=JR, AHTP=8)
						SampleIndex+=1
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IS=IS, IE=IE, SC=SC, JG=JG, IR=IR, JR=JR, AVTP=9, AHTP=9)
						SampleIndex+=1
		for JS in range(1, Ny-11, 4):
			JE=JS+12
			for IG in [0, 1, 2]:
				for IR in [2, 3, 4]:
					for JR in [2, 3, 4]:
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, JS=JS, JE=JE, SC=SC, IG=IG, IR=IR, JR=JR, AVTP=8)
						SampleIndex+=1
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, JS=JS, JE=JE, SC=SC, IG=IG, IR=IR, JR=JR, AHTP=8)
						SampleIndex+=1
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, JS=JS, JE=JE, SC=SC, IG=IG, IR=IR, JR=JR, AVTP=9, AHTP=9)
						SampleIndex+=1
		for IS in range(1, Nx-11, 4):
			IE=IS+12
			for JS in range(1, Ny-7, 4):
				JE=JS+12
				IG=1
				JG=1
				for IR in [2, 3, 4]:
					for JR in [2, 3, 4]:
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AVTP=8)
						SampleIndex+=1
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AHTP=8)
						SampleIndex+=1
						CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AVTP=9, AHTP=9)
						SampleIndex+=1
		for IB in [1, 2, 3, 4]:
			for JG in [0, 1, 2]:
				for IR in [2, 3, 4]:
					for JR in [2, 3, 4]:
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IB=IB, SC=SC, JG=JG, IR=IR, JR=JR, AVTP=8)
						SampleIndex+=1
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IB=IB, SC=SC, JG=JG, IR=IR, JR=JR, AHTP=8)
						SampleIndex+=1
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IB=IB, SC=SC, JG=JG, IR=IR, JR=JR, AVTP=9, AHTP=9)
						SampleIndex+=1
		for JB in [1, 2, 3, 4]:
			for IG in [0, 1, 2]:
				for IR in [2, 3, 4]:
					for JR in [2, 3, 4]:
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, JB=JB, SC=SC, IG=IG, IR=IR, JR=JR, AVTP=8)
						SampleIndex+=1
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, JB=JB, SC=SC, IG=IG, IR=IR, JR=JR, AHTP=8)
						SampleIndex+=1
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, JB=JB, SC=SC, IG=IG, IR=IR, JR=JR, AVTP=9, AHTP=9)
						SampleIndex+=1
	print('Group 2 of Pore-Networks Generated: ', SampleIndex, '---------------------------------')
	# SampleIndex=6474
	# print('Generating 3rd Group of Pore-Networks ------------------------------------------------')
	# for SC in [[3, 2, 1, 0], [5, 6, 7, 8], [0, 1, 5, 6], [8, 7, 3, 2]]:
	# 	for IS in range(1, Nx-15, 4):
	# 		IE=IS+16
	# 		IG=1
	# 		JG=1
	# 		for IR in [1, 2, 3, 4]:
	# 			for JR in [1, 2, 3, 4]:
	# 				CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IS=IS, IE=IE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AVTP=8)
	# 				SampleIndex+=1
	# 				CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IS=IS, IE=IE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AHTP=8)
	# 				SampleIndex+=1
	# 				CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IS=IS, IE=IE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AVTP=9, AHTP=9)
	# 				SampleIndex+=1
	# 	for JS in range(1, Ny-15, 4):
	# 		JE=JS+16
	# 		IG=1
	# 		JG=1
	# 		for IR in [1, 2, 3, 4]:
	# 			for JR in [1, 2, 3, 4]:
	# 				CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, JS=JS, JE=JE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AVTP=8)
	# 				SampleIndex+=1
	# 				CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, JS=JS, JE=JE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AHTP=8)
	# 				SampleIndex+=1
	# 				CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, JS=JS, JE=JE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AVTP=9, AHTP=9)
	# 				SampleIndex+=1
	# 	for IS in range(1, Nx-15, 4):
	# 		IE=IS+16
	# 		for JS in range(1, Ny-15, 4):
	# 			JE=JS+16
	# 			for IG in [0, 1, 2, 3]:
	# 				for IR in [1, 2, 3, 4]:
	# 					for JR in [1, 2, 3, 4]:
	# 						CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, IG=IG, IR=IR, JR=JR, AVTP=8)
	# 						SampleIndex+=1
	# 						CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, IG=IG, IR=IR, JR=JR, AHTP=8)
	# 						SampleIndex+=1
	# 						CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, IG=IG, IR=IR, JR=JR, AVTP=9, AHTP=9)
	# 						SampleIndex+=1
	# 	for IS in range(1, Nx-15, 4):
	# 		IE=IS+16
	# 		for JS in range(1, Ny-15, 4):
	# 			JE=JS+16
	# 			for JG in [0, 1, 2, 3]:
	# 				for IR in [1, 2, 3, 4]:
	# 					for JR in [1, 2, 3, 4]:
	# 						CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, JG=JG, IR=IR, JR=JR, AVTP=8)
	# 						SampleIndex+=1
	# 						CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, JG=JG, IR=IR, JR=JR, AHTP=8)
	# 						SampleIndex+=1
	# 						CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, JG=JG, IR=IR, JR=JR, AVTP=9, AHTP=9)
	# 						SampleIndex+=1
	# 	for B in [1, 2, 3, 4]:
	# 		IB=B
	# 		JB=B
	# 		for RG in [0, 1, 2, 3]:
	# 			for RR in [1, 2, 3, 4]:
	# 				CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IB=IB, JB=JB, SC=SC, RG=RG, RR=RR, AVTP=8)
	# 				SampleIndex+=1
	# 				CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IB=IB, JB=JB, SC=SC, RG=RG, RR=RR, AHTP=8)
	# 				SampleIndex+=1
	# 				CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IB=IB, JB=JB, SC=SC, RG=RG, RR=RR, AVTP=9, AHTP=9)
	# 				SampleIndex+=1
	# print('Group 3 of Pore-Networks Generated: ', SampleIndex, '---------------------------------')
	SampleIndex=9162
	print('Generating 4th Group of Pore-Networks ------------------------------------------------')
	for SC in [[4, 5, 3, 2, 6], [4, 1, 7, 8, 0], [5, 4, 3, 2, 1, 0], [3, 4, 5, 6, 7, 8]]:
		for IS in range(1, Nx-17, 1):
			IE=IS+18
			for JS in range(1, Ny-17, 1):
				JE=JS+18
				for IG in [0, 1]:
					for JG in [0, 1]:
						for IR in [1, 2, 3]:
							for JR in [1, 2, 3]:
								CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AVTP=8)
								SampleIndex+=1
								CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AHTP=8)
								SampleIndex+=1
								CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AVTP=9, AHTP=9)
								# SampleIndex+=1
	print('Group 4 of Pore-Networks Generated: ', SampleIndex, '---------------------------------')

#============================ Main Program ========================================================
# CreatePoreNetworkSamples(Nx= 3, Ny= 3, Folder='/home/xu/work/PoreNetwork0303Samples')
# CreatePoreNetworkSamples(Nx=10, Ny=10, Folder='/home/xu/work/PoreNetwork1010Samples/')
CreatePoreNetworkSamples(Nx=20, Ny=20, Folder='/home/xu/work/PoreNetwork2020Samples/')
# CreatePoreNetworkSamples(Nx=40, Ny=40, Folder='/home/xu/work/PoreNetwork4040Samples/')
# NewNet=pickle.load(open('RandStrNet.net', 'rb'))
