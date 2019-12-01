#==================================================================================================
# Title      : Optimization
# Description: Optimization
# Author     : Zhenyu Xu; westlark@outlook.com
# Start Time : 2019.11.01
# License    : Apache 2.0
#==================================================================================================

#============================ Module Importation ==================================================
import os
import sys
import subprocess
import time

import statistics

import pickle
from copy import deepcopy

from scipy import optimize

# import verb

#==================================================================================================

#============================ Basic Functions =====================================================
# Get the immediate level sub folder names --------------------------------------------------------
def DirectSub(dir):
	return [dir+'/'+name for name in os.listdir(dir) if os.path.isdir(os.path.join(dir, name))]

# Split something to segments according to Splitter
def Segment(Path, Splitter):
	return Path.split(Splitter)

# Create file name --------------------------------------------------------------------------------
def Name(Prefix='Sample', Index=0, Digit=6):
	si=str(Index)
	if Digit<len(si):
		Digit=len(si)
	for i in range(Digit-len(si)):
		si='0'+si
	return Prefix+si



# Get the angle of attack from the last several characters of the folder name ---------------------
def CaseNumber(List, Start=-6):
	p1=List[Start:]
	p2=p1.split('-')[1]
	p3=p2.replace('N', '-')
	p4=float(p3)

	return p4

# Get the case number from the last several characters of the folder name -------------------------
def CaseName(List, Splitter='/', Start=-6):
	p1=List[Start:]
	p2=p1.split(Splitter)[1]

	return p2

# Get the naming and numbering of the given folder ------------------------------------------------
def Separated(List, Splitter='-'):
	SplitIndex=0
	for i in range(len(List)-1, -1, -1):
		if List[i]==Splitter:
			SplitIndex=i
			break
	
	b1=List[:SplitIndex+1]
	p1=List[SplitIndex+1:]
	if Splitter=='-':
		p2=p1.replace('N', '-')
		p3=float(p2)
	else:
		p3=p1

	bp=[]
	bp.append(b1)
	bp.append(p3)

	return bp

# Replace a certain line of a file with given content ---------------------------------------------
def ReadWord(File, LN, WN):
	with open(File,'r') as f:
		line =f.readlines()[LN-1]
		words=line.split()

	f.close()

	return words[WN]

# Replace a certain line of a file with given content ---------------------------------------------
def ReplaceLine(File, LN, Content):
	lines=[]
	with open(File,'r') as f:
		ln=0
		for line in f.readlines():
			ln+=1
			if ln==LN:
				line=Content
			lines.append(line)
	f.close()
	with open(File, 'w') as f:
		for line in lines:
			f.write(line)
	f.close()

	return True

# Check if a process is running or not ------------------------------------------------------------
def Running(PID):
	try:
		os.kill(PID, 0)
	except OSError:
		return False
	else:
		return True

# Check if free CPU available ---------------------------------------------------------------------
def GotFreeCPU():
	CPU=True

	return CPU

# Check the status of RAM usage -------------------------------------------------------------------
def GotFreeRAM():
	RAM=True

	return RAM
#==================================================================================================

#============================ Main Program ========================================================
TotalThread=42
CheckInterval=1

# -------------------------------------------------------------------------------------------------
# WorkDir='/home/xu/work/PoreNetwork2020Samples'
# FixVCount=10424
# RandCount=11396
# WorkDir='/home/xu/work/PoreNetwork1010Samples'
# FixVCount=1138
# RandCount=1462
WorkDir='/home/xu/work/PoreNetwork4040Samples'
FixVCount=1243
RandCount=1266
# WorkDir='/home/xu/work/PoreNetwork1040Samples'
# FixVCount=1117
# RandCount=1225
# WorkDir='/home/xu/work/PoreNetwork4010Samples'
# FixVCount=1117
# RandCount=1225




os.chdir(WorkDir)
SubDirs=DirectSub(WorkDir)
print(SubDirs)

Runs=[]
for subdir in SubDirs:
	FolderName  =Segment(subdir, '/')[-1]
	SettingMatch=Segment(FolderName, '_')[0]
	NetSummary  =Segment(FolderName, '_')[1]
	CasesCount  =FixVCount if NetSummary[1]=='F' else RandCount
	print(FolderName, SettingMatch, NetSummary, CasesCount)
	cmd='cp -f '+SettingMatch+' PoreDy '+subdir
	# subprocess.call(cmd, shell=True)
	print('executed ', cmd)

	cmd='mkdir '+subdir+'/Net'
	# subprocess.call(cmd, shell=True)
	print('executed ', cmd)
	cmd='mv '+subdir+'/*.net '+subdir+'/Net'
	# subprocess.call(cmd, shell=True)
	print('executed ', cmd)

	cmd='mkdir '+subdir+'/EH'
	# subprocess.call(cmd, shell=True)
	print('executed ', cmd)
	cmd='mkdir '+subdir+'/VTK'
	# subprocess.call(cmd, shell=True)
	print('executed ', cmd)
	cmd='mkdir '+subdir+'/Log'
	# subprocess.call(cmd, shell=True)
	print('executed ', cmd)
	cmd='mkdir '+subdir+'/AT'
	# subprocess.call(cmd, shell=True)
	print('executed ', cmd)

	print(subdir)
	os.chdir(subdir)
	for i in range(CasesCount):
	# for i in range(9936,CasesCount, 1):
		jobname=Name(FolderName+'_', i)
		cmd='cp -f '+SettingMatch+' '+jobname
		subprocess.call(cmd, shell=True) # Pay Attention
		print('executed ', cmd)

		cmd='./PoreDy '+jobname+' 1001'
		job=subprocess.Popen(cmd, shell=True) # Pay Attention
		Runs.append(int(job.pid)+1)
		print('executing ', cmd, ' Job ID: ', int(job.pid)+1)
		print(Runs)
		if (i+1)%TotalThread==0:
			time.sleep(1)

		while len(Runs)>=TotalThread:
			for run in Runs:
				if Running(run)==False:
					print('PID ', run, ' Finished!')
					Runs.remove(run)
			time.sleep(CheckInterval)
	
	time.sleep(5)
	cmd='rm '+subdir+'/*.cd'
	subprocess.call(cmd, shell=True) # Pay Attention
	print('executed ', cmd)
	time.sleep(5)
	cmd='mv '+subdir+'/*.at '+subdir+'/AT'
	subprocess.call(cmd, shell=True) # Pay Attention
	print('executed ', cmd)
	time.sleep(5)
	cmd='mv '+subdir+'/*.eh '+subdir+'/EH'
	subprocess.call(cmd, shell=True) # Pay Attention
	print('executed ', cmd)
	time.sleep(5)
	cmd='mv '+subdir+'/*.vtk '+subdir+'/VTK'
	subprocess.call(cmd, shell=True) # Pay Attention
	print('executed ', cmd)
	time.sleep(5)
	cmd='mv '+subdir+'/*.log '+subdir+'/Log'
	subprocess.call(cmd, shell=True) # Pay Attention
	print('executed ', cmd)
	time.sleep(5)

# 	if 
# 	os.chdir(Test)
# 	print('Got in ', Test)
# 	inp=Test+'/cfl3d.inp'
# 	p3d=Test+'/newsurf.p3d'
# 	Freestream=float(ReadWord(File=inp, LN=24, WN=0))
# 	Cases=TopLevelSubDir(Test)
# 	for Case in Cases:
# 		Alpha=Separated(List=Case)[1]
# 		print('   % 7.4f   % 6.3f   % 6.4f % 8.6f   % 6.2f % 7d % 7d' 
# 			% (Freestream, Alpha, 0, 6, 540, 1, 0))
# 		cmd='cp -rf '+inp+' '+p3d+' '+Case+'/'
# 		subprocess.call(cmd, shell=True)
# 		ReplaceLine(File=Case+'/cfl3d.inp', LN=24, 
# 			Content='   % 7.4f   % 6.3f   % 6.4f % 8.6f   % 6.2f % 7d % 7d\n' 
# 			% (Freestream, Alpha, 0, 6, 540, 1, 0))
# 		os.chdir(Case)
# 		cmd='mpirun -n '+str(TotalThread)+' cfl3d_mpi < cfl3d.inp'
# 		print('CFL3D calculating '+Test+'/......-'+str(Alpha)+'......')
# 		subprocess.call(cmd, shell=True)

# pcfl.PostprocessCFL3DProject(WorkDir=WorkDir)


# Study the influence of the ejection speed of airfoil trailing edge ------------------------------
# def StudyEjection(CopyFrom=sys.argv[1], Ejects=[], TotalThread=1):
# 	WorkDir   =Separated(CopyFrom, Splitter='/')[0]
# 	TestFolder=Separated(CopyFrom, Splitter='/')[1]
# 	Name      =Separated(TestFolder, Splitter='-')[0]
# 	Number    =Separated(TestFolder, Splitter='-')[1]

# 	Runs=[]

# 	for Ejection in Ejects:
# 		CopyTo=WorkDir+Name+str(Ejection)
# 		cmd='cp -rf '+CopyFrom+' '+CopyTo
# 		subprocess.call(cmd, shell=True)
# 		Cases=TopLevelSubDir(CopyTo)
# 		for Case in Cases:
# 			os.chdir(Case)
# 			print(Case)

# 			ReplaceLine(File=Case+'/cfl3d.inp', LN=2, Content='cfl3d.xyz\n')
# 			RatioT=1+(1.4-1)/2*pow(Ejection, 2)
# 			RatioP=pow(RatioT, 1.4/(1.4-1))
# 			ReplaceLine(File=Case+'/cfl3d.inp', LN=248, Content=
# 				"\t%4.2f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n" % (Ejection, RatioP, RatioT, 0, 0))
# 			job=subprocess.Popen('cfl3d_seq<cfl3d.inp &', shell=True)
# 			Runs.append(int(job.pid)+1)
# 			print(int(job.pid)+1)
# 			print(Runs)
# 			time.sleep(10)

# 			while len(Runs)>=TotalThread:
# 				for run in Runs:
# 					if Running(run)==False:
# 						print('kill ', run)
# 						Runs.remove(run)
# 				time.sleep(CheckInterval)

# 	return True