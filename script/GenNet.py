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
MyDPI=144
OpenMap ={'S':['S'], 'D':['S', 'N'], 'C':['E', 'N'], 'T':['E', 'N', 'W'], 'Q':['E', 'N', 'W', 'S']}
CrossMap={False:'S', True:'C'}
FixVMap ={False:'R', True:'F'}
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
			self.Pool       [ThroatType[i]]=ThroatCount[i]
			self.InitialPool[ThroatType[i]]=ThroatCount[i]
			self.InitialTotal             +=ThroatCount[i]

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
	def FetchOne(self, Index=0, Sub=[], Exc=[], Method='Random', Install=1):
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
			return Selection*Install

	# Return One ----------------------------------------------------------------------------------
	def ReturnOne(self, T):
		T=abs(T)
		if self.InPool(T):
			self.Pool[T]+=1
			return True
		else:
			print('Warning: -------------------------------cannot return, no throat type matched!')
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
		if 'W' in self.Open:
			self.VTRange[0][0] =1
			self.StrTCount    -=self.Ny
		if 'E' in self.Open:
			self.VTRange[0][1] =self.Nx
			self.StrTCount    -=self.Ny
		if 'S' in self.Open:
			self.HTRange[1][0] =1
			self.StrTCount    -=self.Nx
		if 'N' in self.Open:
			self.HTRange[1][1] =self.Ny
			self.StrTCount    -=self.Nx
		self.TotTCount=self.StrTCount+self.CrsTCount

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
			print('Warning: -------------------------Straight Throat Pool Created but Not Enough!')
		if self.Cross and self.CrsTCount>self.CrsTPool.ExistTotal():
			print('Warning: ----------------------------Cross Throat Pool Created but Not Enough!')

	# set a new name ------------------------------------------------------------------------------
	def SetName(self, Name):
		self.Name=Name
		return True

	# Output the file binary and other related files ----------------------------------------------
	def Output(self, Folder='', Disp=False, Option=['Dump', 'Write', 'Gray', 'RGB']):
		if 'Dump' in Option: pickle.dump(self, open(Folder+self.Name+'.net', 'wb'))
		if Disp: print('Successfully saved pore-network binary file.')

		if 'Write' in Option or 'Gray' in Option or 'RGB' in Option:
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

			if 'Write' in Option:
				with open(Folder+self.Name+'.at', 'w') as wat:
					for j in range(My-1, -1, -1):
						for i in range(Mx):
							if Matrix[j][i]==0:
								wat.write('% 9.6e\t% 9d\t' % (0, 0))
							else:
								throattype=Matrix[j][i]
								flip=throattype/abs(throattype)
								row =abs(throattype)%10
								col =abs(throattype)//10
								wat.write('% 9.6e\t% 9d\t' % (ThroatChoice[row][col].D*flip, ThroatChoice[row][col].S))
						wat.write('\n')

			if 'Gray' in Option:
				Dia=np.zeros((My, Mx), dtype=float)
				for j in range(My):
					for i in range(Mx):
						if Matrix[j][i]==0:
							pass
						else:
							throattype=Matrix[j][i]
							flip      =throattype/abs(throattype)
							row       =abs(throattype)%10
							col       =abs(throattype)//10
							Dia[j][i] =DToGray[ThroatChoice[row][col].D*flip]
				plt.imshow(Dia, cmap="gray", vmin=min(DToGray.values()), vmax=max(DToGray.values()))
				if Disp: plt.show()
				plt.savefig(Folder+self.Name+'.png')

			if 'RGB' in Option:
				pass

		return True

	# Get Throat info at a certain position -------------------------------------------------------
	def GetT(self, TP, I, J):
		if TP=='VT':
			if (I>=self.VTRange[0][0] and I<self.VTRange[0][1] and J>=self.VTRange[1][0] and J<self.VTRange[1][1]):
				if self.StrTPool.InPool(self.VT[J][I]):
					return ['Exist&Assigned', self.VT[J][I]]
				else:
					return ['Exist&NotAssigned']
			else:
				return ['NotExist']
		elif TP=='HT':
			if (I>=self.HTRange[0][0] and I<self.HTRange[0][1] and J>=self.HTRange[1][0] and J<self.HTRange[1][1]):
				if self.StrTPool.InPool(self.HT[J][I]):
					return ['Exist&Assigned', self.HT[J][I]]
				else:
					return ['Exist&NotAssigned']
			else:
				return ['NotExist']
		elif TP=='CT':
			if (I>=self.CTRange[0][0] and I<self.CTRange[0][1] and J>=self.CTRange[1][0] and J<self.CTRange[1][1]) and self.Cross:
				if self.CrsTPool.InPool(self.CT[J][I]):
					return ['Exist&Assigned', self.CT[J][I]]
				else:
					return ['Exist&NotAssigned']
			else:
				return ['NotExist']
		else:
			return ['TotallyWrong']

	# Return Throat at a certain position back to appropriate pool --------------------------------
	def ReturnT(self, TP, I, J):
		InfoT=self.GetT(TP, I, J)
		if InfoT[0]=='Exist&Assigned':
			if   TP=='VT':
				self.StrTPool.ReturnOne(InfoT[1])
				self.VT[J][I]=0
			elif TP=='HT':
				self.StrTPool.ReturnOne(InfoT[1])
				self.HT[J][I]=0
			elif TP=='CT':
				self.CrsTPool.ReturnOne(InfoT[1])
				self.CT[J][I]=0
			else:
				print('Warning: No Throat to Return!')
				return False
			return True
		else:
			return False

	# Assign Throat at a certain position ---------------------------------------------------------
	def AssignT(self, TP, I, J, Index=0, Sub=[], Exc=[], Method='Random', Install=1):
		if   TP=='VT' and I>=self.VTRange[0][0] and I<self.VTRange[0][1] and J>=self.VTRange[1][0] and J<self.VTRange[1][1]:
			self.ReturnT(TP, I, J)
			self.VT[J][I]=self.StrTPool.FetchOne(Index, Sub, Exc, Method, Install)
		elif TP=='HT' and I>=self.HTRange[0][0] and I<self.HTRange[0][1] and J>=self.HTRange[1][0] and J<self.HTRange[1][1]:
			self.ReturnT(TP, I, J)
			self.HT[J][I]=self.StrTPool.FetchOne(Index, Sub, Exc, Method, Install)
		elif TP=='CT' and I>=self.CTRange[0][0] and I<self.CTRange[0][1] and J>=self.CTRange[1][0] and J<self.CTRange[1][1] and self.Cross:
			self.ReturnT(TP, I, J)
			self.CT[J][I]=self.CrsTPool.FetchOne(Index, Sub, Exc, Method, Install)
		else:
			print('Warning: -------------------------------------------------AssignT not succeed!', TP, I, J, self.VTRange, self.HTRange, self.CTRange)
			return False
		return True

	# Assign Throat distribution in a box of certain position -------------------------------------
	def AssignBox(self, TP='',
		                Start=[0, 0], End=[0, 0], Band=[0, 0], , Jump=[0, 0]
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
		if Start[0]<Range[0][0]: Start[0]=Range[0][0]
		if Start[1]<Range[1][0]: Start[1]=Range[1][0]
		if End  [0]>Range[0][1]: End  [0]=Range[0][1]
		if End  [1]>Range[1][1]: End  [1]=Range[1][1]
		if Sub==[]:
			Sub=[key for key in PoolChoice.Pool.keys()]
			Exc=[]
		else:
			Exc=Sub
		# print('AssignBox Before: ', TP, Start, End, Sub, Exc)
		StrAdd       =0
		StrReplace   =0
		CrsAdd       =0
		CrsReplace   =0
		StrPoolBefore=self.StrTPool.ExistTotal()
		if self.Cross: CrsPoolBefore=self.CrsTPool.ExistTotal()

		for I in range(Start[0], End[0], 1+Jump[0]):
			for J in range(Start[1], End[1], 1+Jump[1]):
				if Band==[0, 0] or (\
				   (Band[0]> 0 and Band[1]> 0) and not(I>=Start[0]+Band[0] and I<End[0]-Band[0] and \
				                                       J>=Start[1]+Band[1] and J<End[1]-Band[1])) or (\
				   (Band[0]> 0 and Band[1]==0) and not(I>=Start[0]+Band[0] and I<End[0]-Band[0])) or (\
				   (Band[0]==0 and Band[1]> 0) and not(J>=Start[1]+Band[1] and J<End[1]-Band[1])):

					InfoT=self.GetT(TP, I, J)
					if InfoT[0]=='Exist&Assigned':
						if TP=='VT' or TP=='HT':
							StrReplace+=1
						elif TP=='CT':
							CrsReplace+=1
					elif InfoT[0]=='Exist&NotAssigned':
						if TP=='VT' or TP=='HT':
							StrAdd    +=1
						elif TP=='CT':
							CrsAdd    +=1

					for ring in range(min(Band[0], Band[1])+1):
						if not (I>=Start[0]+ring+1 and I<End[0]-ring-1 and \
						        J>=Start[1]+ring+1 and J<End[1]-ring-1) \
						   and (I>=Start[0]+ring   and I<End[0]-ring and \
						        J>=Start[1]+ring   and J<End[1]-ring):
							RingIndex=ring
					Selection=Pick(List=Sub, Index=Grad[0]*(I-Start[0])//Repeat[0]+\
					                               Grad[1]*(J-Start[1])//Repeat[1]+\
					                               Grad[2]* RingIndex  //Repeat[2], Method='Rotate')
					Install=1
					if TP=='VT' or TP=='HT':
						pass
					elif TP=='CT':
						if Flip==0:
							Install=Pick(List=[-1, 1])
						elif Flip==2 and I%2==0:
							Install=-1
						elif Flip==3 and J%2==0:
							Install=-1
						elif Flip==-1 or Flip==1: # Flip= -1 or 1
							Install=Flip
					T=self.AssignT(TP, I, J, Selection, Sub, Exc, 'Rotate', Install)
		if StrAdd!=StrPoolBefore-self.StrTPool.ExistTotal():
			print('AssignBox Error: Straight Throat added in Net Not Equal diminished in Pool!')
		if self.Cross:
			if CrsAdd!=CrsPoolBefore-self.CrsTPool.ExistTotal():
				print('AssignBox Error: Cross    Throat added in Net Not Equal diminished in Pool!')
		return StrAdd+StrReplace+CrsAdd+CrsReplace

	# Slice a Region from original Matrix ---------------------------------------------------------
	def RandomRest(self):
		vt=0
		ht=0
		StrPoolBefore=self.StrTPool.ExistTotal()
		ct=0
		if self.Cross: CrsPoolBefore=self.CrsTPool.ExistTotal()
		for I in range(self.VTRange[0][0], self.VTRange[0][1], 1):
			for J in range(self.VTRange[1][0], self.VTRange[1][1], 1):
				InfoT=self.GetT('VT', I, J)
				if InfoT[0]=='Exist&NotAssigned':
					self.AssignT(TP='VT', I=I, J=J)
					vt+=1
		for I in range(self.HTRange[0][0], self.HTRange[0][1], 1):
			for J in range(self.HTRange[1][0], self.HTRange[1][1], 1):
				InfoT=self.GetT('HT', I, J)
				if InfoT[0]=='Exist&NotAssigned':
					self.AssignT(TP='HT', I=I, J=J)
					ht+=1
		if self.Cross:
			for I in range(self.CTRange[0][0], self.CTRange[0][1], 1):
				for J in range(self.CTRange[1][0], self.CTRange[1][1], 1):
					InfoT=self.GetT('CT', I, J)
					if InfoT[0]=='Exist&NotAssigned':
						self.AssignT(TP='CT', I=I, J=J)
						ct+=1
		if vt+ht!=StrPoolBefore-self.StrTPool.ExistTotal():
			print('RandomRest Str Error: ', vt+ht+ct, '\tVT: ', vt, '\tHT: ', ht, '\tCT: ', ct)
		
		if self.Cross:
			if ct   !=CrsPoolBefore-self.CrsTPool.ExistTotal():
				print('RandomRest Str Error: ', vt+ht+ct, '\tVT: ', vt, '\tHT: ', ht, '\tCT: ', ct)
		return vt+ht+ct

	# Check if it is correct ----------------------------------------------------------------------
	def Check(self):
		vt=0
		ht=0
		ct=0
		for I in range(self.Nx+1):
			for J in range(self.Ny+1):
				for TP in ['VT', 'HT', 'CT']:
					InfoT=self.GetT(TP=TP, I=I, J=J)
					if   TP=='VT':
						if I>=self.VTRange[0][0] and I<self.VTRange[0][1] and J>=self.VTRange[1][0] and J<self.VTRange[1][1]:
							if InfoT[0]=='Exist&Assigned':
								vt+=1
							else:
								print('Error: ------------------In Range Not Assigned: ', TP, I, J)
						else:
							if InfoT[0]=='NotExist' and self.VT[J][I]==0:
								pass
							else:
								print('Error: --------------Not In Range But Assigned: ', TP, I, J)

					elif TP=='HT':
						if I>=self.HTRange[0][0] and I<self.HTRange[0][1] and J>=self.HTRange[1][0] and J<self.HTRange[1][1]:
							if InfoT[0]=='Exist&Assigned':
								ht+=1
							else:
								print('Error: ------------------In Range Not Assigned: ', TP, I, J)
						else:
							if InfoT[0]=='NotExist' and self.HT[J][I]==0:
								pass
							else:
								print('Error: --------------Not In Range But Assigned: ', TP, I, J)
					elif TP=='CT':
						if I>=self.CTRange[0][0] and I<self.CTRange[0][1] and J>=self.CTRange[1][0] and J<self.CTRange[1][1] and self.Cross:
							if InfoT[0]=='Exist&Assigned':
								ct+=1
							else:
								print('Error: ------------------In Range Not Assigned: ', TP, I, J)
						else:
							if InfoT[0]=='NotExist' and self.CT[J][I]==0:
								pass
							else:
								print('Error: --------------Not In Range But Assigned: ', TP, I, J)
		if vt+ht!=self.StrTPool.InitialTotal-self.StrTPool.ExistTotal():
			print('Error: -------Total Straight Throat in Net is not equal to diminished in Pool!')
		if self.Cross:
			if ct   !=self.CrsTPool.InitialTotal-self.CrsTPool.ExistTotal():
				print('Error: -------Total Cross    Throat in Net is not equal to diminished in Pool!')
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
# Create samples for One Assignment Region for each of the 4 kinds of networks ====================
def Create(Nx=20, Ny=20, Folder='', SIN=0, OutOpt=['Dump', 'Write'], \
	       Cross=[False, True], FixV=[False, True], AVTP=0, AHTP=0, ACTP=0, \
	       VIS=0, VJS=0, VIE=0, VJE=0, VIB=0, VJB=0, VIJ=0, VJJ=0, SplitRangeV=False, \
	       HIS=0, HJS=0, HIE=0, HJE=0, HIB=0, HJB=0, HIJ=0, HJJ=0, SplitRangeH=False, \
	       CIS=0, CJS=0, CIE=0, CJE=0, CIB=0, CJB=0, CIJ=0, CJJ=0, SplitRangeC=False, \
	       SubV=[], SubH=[], SubC=[], \
	       VIG=0, VJG=0, VRG=0, HIG=0, HJG=0, HRG=0, CIG=0, CJG=0, CRG=0, \
	       VIR=1, VJR=1, VRR=1, HIR=1, HJR=1, HRR=1, CIR=1, CJR=1, CRR=1, \
	       CFlip=0):
	for openkey, openvalue in OpenMap.items():
		for cross in Cross:
			for fixv in FixV:
				# ---------------------------------------------------------------------------------
				prefix='PN_'+CrossMap[cross]+FixVMap[fixv]+openkey
				name  =Name(Prefix=prefix+'_', Index=SIN)
				print('Creating ', name, ':')
				PN=PoreNetwork(Name=name, \
				               Nx=Nx, Ny=Ny, Open=openvalue, Cross=cross , FixV=fixv, \
				               StrTType=TT, CrsTType=TT)
				if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<AVTP:
					STAND=[[[VIS, VJS], [VIE, VJE]]]
					if SplitRangeV:
						STAND=Pick(List=[[[[VIS, 0], [VIE, 0]], [[0, VJS], [0, VJE]]], [[[0, VJS], [0, VJE]], [[VIS, 0], [VIE, 0]]]])
					for stnd in STAND:
						PN.AssignBox(TP='VT', Start=stnd[0], End=stnd[1], Band=[VIB, VJB], Jump=[VIJ, VJJ]\
						             Sub=SubV, Grad=[VIG, VJG, VRG], Repeat=[VIR, VJR, VRR])
						print('\t', stnd[0], stnd[1], [VIB, VJB], SubV, [VIG, VJG, VRG], [VIR, VJR, VRR])
				if Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<AHTP:
					STAND=[[[HIS, HJS], [HIE, HJE]]]
					if SplitRangeH:
						STAND=Pick(List=[[[[HIS, 0], [HIE, 0]], [[0, HJS], [0, HJE]]], [[[0, HJS], [0, HJE]], [[HIS, 0], [HIE, 0]]]])
					for stnd in STAND:
						PN.AssignBox(TP='HT', Start=stnd[0], End=stnd[1], Band=[HIB, HJB], Jump=[HIJ, HJJ]\
						             Sub=SubH, Grad=[HIG, HJG, HRG], Repeat=[HIR, HJR, HRR])
						print('\t', stnd[0], stnd[1], [HIB, HJB], SubH, [HIG, HJG, HRG], [HIR, HJR, HRR])
				if cross and Pick(List=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])<ACTP:
					STAND=[[[CIS, CJS], [CIE, CJE]]]
					if SplitRangeC:
						STAND=Pick(List=[[[[CIS, 0], [CIE, 0]], [[0, CJS], [0, CJE]]], [[[0, CJS], [0, CJE]], [[CIS, 0], [CIE, 0]]]])
					for stnd in STAND:
						PN.AssignBox(TP='CT', Start=stnd[0], End=stnd[1], Band=[CIB, CJB], Jump=[CIJ, CJJ]\
						             Sub=SubC, Grad=[CIG, CJG, CRG], Repeat=[CIR, CJR, CRR], \
						             Flip=CFlip)
						print('\t', stnd[0], stnd[1], [CIB, CJB], SubC, [CIG, CJG, CRG], [CIR, CJR, CRR])
				PN.RandomRest()
				PN.Output(Folder+prefix+'/', Option=OutOpt)
				PN.Check()
				print('--------------------------------------------------------------------------')
				# ---------------------------------------------------------------------------------
	return True
# =================================================================================================
def CreatePoreNetworkSamples(Nx=20, Ny=20, Folder='', OutOpt=['Dump', 'Write']):
	print('Prepare the folders ------------------------------------------------------------------')
	cmd='mkdir '+Folder
	subprocess.call(cmd, shell=True)
	for openkey, openvalue in OpenMap.items():
		for cross in Cross:
			for fixv in FixV:
				cmd='mkdir '+Folder+CrossMap[cross]+FixVMap[fixv]+openkey+'/'
				subprocess.call(cmd, shell=True)

	print('Initialize Sample Index --------------------------------------------------------------')
	SampleIndex    =0
	RandBreakPoints=[]
	RandSampleIndex=SampleIndex
	RandBreakPoints.append(RandSampleIndex)
	FixVBreakPoints=[]
	FixVSampleIndex=SampleIndex
	FixVBreakPoints.append(FixVSampleIndex)

	print('Generating Pore-Network Samples ------------------------------------------------------')
	for sub1 in [[15], [25], [35], [45], [55], [65], [75], [85], [95]]:
		for sub2 in [[15], [25], [35], [45], [55], [65], [75], [85], [95]]:
			for jump in [0, 1, 2, 3, 4]:
				Create(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=RandSampleIndex, FixV=[False], AVTP=10, VIJ=jump, SubV=sub1)
				RandSampleIndex+=1
				Create(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=RandSampleIndex, FixV=[False], AHTP=10, HJJ=jump, SubH=sub1)
				RandSampleIndex+=1
				Create(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=RandSampleIndex, FixV=[False], Cross=[False], AVTP=10, SubV=sub1, AHTP=10, HJJ=jump, SubC=sub2)
				Create(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=RandSampleIndex, FixV=[False], Cross=[True ], AVTP=10, SubV=sub1, ACTP=10, CJJ=jump, SubC=sub2, CFlip=2)
				RandSampleIndex+=1
				Create(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=RandSampleIndex, FixV=[False], Cross=[False], AHTP=10, SubV=sub1, AVTP=10, VJJ=jump, SubC=sub2)
				Create(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=RandSampleIndex, FixV=[False], Cross=[True ], AHTP=10, SubH=sub1, ACTP=10, CIJ=jump, SubC=sub2, CFlip=3)
				RandSampleIndex+=1
	RandBreakPoints.append(RandSampleIndex)
	print('Pore-Network Samples Random Group No.', len(RandBreakPoints)-1, ': ', RandBreakPoints[-1], '---------------------')
	# FixVBreakPoints.append(FixVSampleIndex)
	# print('Pore-Network Samples Fixed  Group No.', len(FixVBreakPoints)-1, ': ', FixVBreakPoints[-1], '---------------------')

	for b in range(1, min(Nx, Ny)//2, 1):

	for sub1 in [[15, 25], [35, 45], [75, 65], [95, 85], [25, 55], [45, 55], [65, 55], [85, 55]]:
		for sub2 in [[15, 25, 35, 45, 55, 65, 75, 85, 95], [95, 85, 75, 65, 55, 45, 35, 25, 15]]:
			for r1 in [1, 2, 3, 4]:
				for r2 in [1, 2, 3, 4]:
					for IS in range(0, Nx-7, 4):
						IE=IS+8

						Create(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, AVTP=10, SubV=sub1, VIG=1, VIR=r1)
						SampleIndex+=1
						Create(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, AVTP=10, SubV=sub1, VIG=1, VIR=r1)
						SampleIndex+=1
						Create(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IS=IS, IE=IE, SC=SC, IG=IG, IR=IR, JR=JR, AVTP=9, AHTP=9)
						SampleIndex+=1
			for JS in range(0, Ny-7, 4):
				JE=JS+8
				for JG in [0, 1]:
					for IR in [2, 3]:
						for JR in [1, 3]:
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, JS=JS, JE=JE, SC=SC, JG=JG, IR=IR, JR=JR, AVTP=8)
						SampleIndex+=1
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, JS=JS, JE=JE, SC=SC, JG=JG, IR=IR, JR=JR, AHTP=8)
						SampleIndex+=1
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, JS=JS, JE=JE, SC=SC, JG=JG, IR=IR, JR=JR, AVTP=9, AHTP=9)
						SampleIndex+=1
		for IS in range(1, Nx-3, 4):
			IE=IS+4
			for JS in range(1, Ny-3, 4):
				JE=JS+4
				IG=1
				JG=1
				for IR in [1, 2, 3, 4]:
					for JR in [1, 2, 3, 4]:
						CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AVTP=8)
						SampleIndex+=1
						CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AHTP=8)
						SampleIndex+=1
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AVTP=9, AHTP=9)
						SampleIndex+=1
		for IB in [1, 2, 3, 4]:
			for IG in [0, 1]:
				for IR in [2, 4]:
					for JR in [1, 3]:
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IB=IB, SC=SC, IG=IG, IR=IR, JR=JR, AVTP=8)
						SampleIndex+=1
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IB=IB, SC=SC, IG=IG, IR=IR, JR=JR, AHTP=8)
						SampleIndex+=1
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IB=IB, SC=SC, IG=IG, IR=IR, JR=JR, AVTP=9, AHTP=9)
						SampleIndex+=1
		for JB in [1, 2, 3, 4]:
			for JG in [0, 1]:
				for IR in [1, 3]:
					for JR in [2, 4]:
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, JB=JB, SC=SC, JG=JG, IR=IR, JR=JR, AVTP=8)
						SampleIndex+=1
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, JB=JB, SC=SC, JG=JG, IR=IR, JR=JR, AHTP=8)
						SampleIndex+=1
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, JB=JB, SC=SC, JG=JG, IR=IR, JR=JR, AVTP=9, AHTP=9)
						SampleIndex+=1
	print('Group 1 of Pore-Networks Generated: ', SampleIndex, '---------------------------------')
	SampleIndex=3072
	print('Generating 2nd Group of Pore-Networks ------------------------------------------------')
	for SC in [[15, 25, 35], [55, 45, 65], [95, 85, 75]]:
		for IS in range(1, Nx-11, 4):
			IE=IS+12
			for JG in [0, 1, 2]:
				for IR in [2, 3, 4]:
					for JR in [2, 3, 4]:
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IS=IS, IE=IE, SC=SC, JG=JG, IR=IR, JR=JR, AVTP=8)
						SampleIndex+=1
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IS=IS, IE=IE, SC=SC, JG=JG, IR=IR, JR=JR, AHTP=8)
						SampleIndex+=1
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IS=IS, IE=IE, SC=SC, JG=JG, IR=IR, JR=JR, AVTP=9, AHTP=9)
						SampleIndex+=1
		print('Generated 2.1 Group of Pore-Networks ------------------------------------------------', SampleIndex)
		for JS in range(1, Ny-11, 4):
			JE=JS+12
			for IG in [0, 1, 2]:
				for IR in [2, 3, 4]:
					for JR in [2, 3, 4]:
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, JS=JS, JE=JE, SC=SC, IG=IG, IR=IR, JR=JR, AVTP=8)
						SampleIndex+=1
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, JS=JS, JE=JE, SC=SC, IG=IG, IR=IR, JR=JR, AHTP=8)
						SampleIndex+=1
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, JS=JS, JE=JE, SC=SC, IG=IG, IR=IR, JR=JR, AVTP=9, AHTP=9)
						SampleIndex+=1
		print('Generated 2.2 Group of Pore-Networks ------------------------------------------------', SampleIndex)
		for IS in range(1, Nx-11, 4):
			IE=IS+12
			for JS in range(1, Ny-11, 4):
				JE=JS+12
				IG=1
				JG=1
				for IR in [2, 3, 4]:
					for JR in [2, 3, 4]:
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AVTP=8)
						SampleIndex+=1
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AHTP=8)
						SampleIndex+=1
						CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AVTP=9, AHTP=9)
						SampleIndex+=1
		print('Generated 2.3 Group of Pore-Networks ------------------------------------------------', SampleIndex)
		for IB in [1, 2, 3, 4]:
			for JG in [0, 1, 2]:
				for IR in [2, 3, 4]:
					for JR in [2, 3, 4]:
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IB=IB, SC=SC, JG=JG, IR=IR, JR=JR, AVTP=8)
						SampleIndex+=1
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IB=IB, SC=SC, JG=JG, IR=IR, JR=JR, AHTP=8)
						SampleIndex+=1
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IB=IB, SC=SC, JG=JG, IR=IR, JR=JR, AVTP=9, AHTP=9)
						SampleIndex+=1
		print('Generated 2.4 Group of Pore-Networks ------------------------------------------------', SampleIndex)
		for JB in [1, 2, 3, 4]:
			for IG in [0, 1, 2]:
				for IR in [2, 3, 4]:
					for JR in [2, 3, 4]:
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, JB=JB, SC=SC, IG=IG, IR=IR, JR=JR, AVTP=8)
						SampleIndex+=1
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, JB=JB, SC=SC, IG=IG, IR=IR, JR=JR, AHTP=8)
						SampleIndex+=1
						CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, JB=JB, SC=SC, IG=IG, IR=IR, JR=JR, AVTP=9, AHTP=9)
						SampleIndex+=1
	print('Group 2 of Pore-Networks Generated: ', SampleIndex, '---------------------------------')
	SampleIndex=6474
	print('Generating 3rd Group of Pore-Networks ------------------------------------------------')
	for SC in [[45, 35, 25, 15], [65, 75, 85, 95], [15, 25, 65, 75], [95, 85, 45, 35]]:
		for IS in range(1, Nx-15, 4):
			IE=IS+16
			IG=1
			JG=1
			for IR in [1, 2, 3, 4]:
				for JR in [1, 2, 3, 4]:
					CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IS=IS, IE=IE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AVTP=8)
					SampleIndex+=1
					CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IS=IS, IE=IE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AHTP=8)
					SampleIndex+=1
					CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IS=IS, IE=IE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AVTP=9, AHTP=9)
					SampleIndex+=1
		for JS in range(1, Ny-15, 4):
			JE=JS+16
			IG=1
			JG=1
			for IR in [1, 2, 3, 4]:
				for JR in [1, 2, 3, 4]:
					CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, JS=JS, JE=JE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AVTP=8)
					SampleIndex+=1
					CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, JS=JS, JE=JE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AHTP=8)
					SampleIndex+=1
					CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, JS=JS, JE=JE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AVTP=9, AHTP=9)
					SampleIndex+=1
		for IS in range(1, Nx-15, 4):
			IE=IS+16
			for JS in range(1, Ny-15, 4):
				JE=JS+16
				for IG in [0, 1, 2, 3]:
					for IR in [1, 2, 3, 4]:
						for JR in [1, 2, 3, 4]:
							CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, IG=IG, IR=IR, JR=JR, AVTP=8)
							SampleIndex+=1
							CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, IG=IG, IR=IR, JR=JR, AHTP=8)
							SampleIndex+=1
							CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, IG=IG, IR=IR, JR=JR, AVTP=9, AHTP=9)
							SampleIndex+=1
		for IS in range(1, Nx-15, 4):
			IE=IS+16
			for JS in range(1, Ny-15, 4):
				JE=JS+16
				for JG in [0, 1, 2, 3]:
					for IR in [1, 2, 3, 4]:
						for JR in [1, 2, 3, 4]:
							CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, JG=JG, IR=IR, JR=JR, AVTP=8)
							SampleIndex+=1
							CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, JG=JG, IR=IR, JR=JR, AHTP=8)
							SampleIndex+=1
							CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, JG=JG, IR=IR, JR=JR, AVTP=9, AHTP=9)
							SampleIndex+=1
		for B in [1, 2, 3, 4]:
			IB=B
			JB=B
			for RG in [0, 1, 2, 3]:
				for RR in [1, 2, 3, 4]:
					CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IB=IB, JB=JB, SC=SC, RG=RG, RR=RR, AVTP=8)
					SampleIndex+=1
					CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IB=IB, JB=JB, SC=SC, RG=RG, RR=RR, AHTP=8)
					SampleIndex+=1
					CreateOne(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IB=IB, JB=JB, SC=SC, RG=RG, RR=RR, AVTP=9, AHTP=9)
					SampleIndex+=1
	print('Group 3 of Pore-Networks Generated: ', SampleIndex, '---------------------------------')
	SampleIndex=9162
	print('Generating 4th Group of Pore-Networks ------------------------------------------------')
	for SC in [[55, 65, 45, 35, 75], [55, 25, 85, 95, 15], [65, 55, 45, 35, 25, 15], [45, 55, 65, 75, 85, 95]]:
		for IS in range(1, Nx-17, 1):
			IE=IS+18
			for JS in range(1, Ny-17, 1):
				JE=JS+18
				for IG in [0, 1]:
					for JG in [0, 1]:
						for IR in [1, 2, 3]:
							for JR in [1, 2, 3]:
								CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AVTP=8)
								SampleIndex+=1
								CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AHTP=8)
								SampleIndex+=1
								CreateTwo(Nx=Nx, Ny=Ny, Folder=Folder, OutOpt=OutOpt, SIN=SampleIndex, IS=IS, IE=IE, JS=JS, JE=JE, SC=SC, IG=IG, JG=JG, IR=IR, JR=JR, AVTP=9, AHTP=9)
								# SampleIndex+=1
	print('Group 4 of Pore-Networks Generated: ', SampleIndex, '---------------------------------')
	print('Generating 5th Group of Pore-Networks ------------------------------------------------')
	# for SC in [[15], [25], [35], [45], [55], [65], [75], [85], [95]]:


#============================ Main Program ========================================================
# CreatePoreNetworkSamples(Nx= 3, Ny= 3, Folder='/home/xu/work/PoreNetwork0303Samples')
# CreatePoreNetworkSamples(Nx=10, Ny=10, Folder='/home/xu/work/PoreNetwork1010Samples/')
CreatePoreNetworkSamples(Nx=20, Ny=20, Folder='/home/xu/work/PoreNetwork2020Samples/')
# CreatePoreNetworkSamples(Nx=40, Ny=40, Folder='/home/xu/work/PoreNetwork4040Samples/')
# NewNet=pickle.load(open('RandStrNet.net', 'rb'))
