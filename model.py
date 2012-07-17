#!/usr/bin/python

import sys
import os
import math
import time
import random
import settings
import copy
import numpy as np
from utils.tools import *

class seqDNA:
        length = 0
        seq = []
        site = []
	siteNum = 0
        name = ""
        def __init__(self):
                pass

        def __init__(self, name, seq=[]):
                self.name = name
                self.insertSeq(seq)
	''' 
	   @parameter
		seq: string or list; the fragment need to insert into the seq
		site: int or list;the site index of the fragment, if not given, the site will be incremeted from the last site
		start: int;the insert position in the previous sequence.If not given, when the sequence is empty, the function 
			   used to assign the fragment to the sequence. While the sequence is unempty, the fragment is appended
			    into the sequence.
	   @return
		when the type of parameter seq is wrong, the program is halted.
	'''
        def insertSeq(self, seq, site=[], start=-1):
                if self.length < 1:
                        start = 0 
                elif -1 == start:
                        start = self.length
                if isinstance(seq, str):
                        seq = list(seq)
                else: 
                        assert type(seq) is list, "The parameter seq for insertSeq wrong %s."
                newLen = len(seq)
                newSite = self.siteNum + 1 
		self.siteNum += newLen
		self.length += newLen
                self.seq = self.seq[:start] + seq + self.seq[start:]
		if isinstance(site, int):
			site=[site]
		elif len(site) < 1:
			site = range(newSite, self.siteNum+1)
                self.site = self.site[:start] + site + self.site[start:]

	def removeGap(self):
		delArr = []
		for i in range(self.length):
			if settings.GAP == self.seq[i]:
				delArr.append(i)
		delNum = 0
		for i in delArr:
			del self.seq[i - delNum] 
			del self.site[i -delNum]
			delNum += 1
				
		self.length -= len(delArr)

class eModel:
	@staticmethod
	def emission(t1, t2, a=None):
		pass

	@staticmethod
	def emissionProb(a,b,t1, t2):
		return 1.0

		
class dnaModelJC69(eModel):

	modelName = "JC69"	
	# remain probability	
	@staticmethod
	def p0(t):
		return (1.0/4 + 3.0/4 * math.exp(-4 * settings.LAMBDA * (t/3))) 

	# substitution probability
	@staticmethod
	def p1(t):
		return (1.0/4 - 1.0/4 * math.exp(-4 * settings.LAMBDA * (t/3)))

	eMatrix = None	
	
	@staticmethod
	def setEmissionProb(t):
		dnaModelJC69.eMatrix = [[dnaModelJC69.p1(t) for i in range(4)] for j in range(4)]
		for i in range(4):
			dnaModelJC69.eMatrix[i][i] = dnaModelJC69.p0(t)


	@staticmethod
	def emission(t1, t2, a=None):
		if None == a:
			a = sampleNucleotide()	
	#	pBaseSub = [dnaModelJC69.p1(t1+t2) for i in range(4)]
	#	pBaseSub[settings.BASE.index(a)] = dnaModelJC69.p0(t1+t2) 	
		y = settings.BASE[randomIndex(dnaModelJC69.eMatrix[settings.BASE.index(a)])]
		'''
		pBaseSub = [dnaModelJC69.p1(t2) for i in range(4)]
		pBaseSub[settings.BASE.index(a)] = dnaModelJC69.p0(t2) 	
		a = y
		y = settings.BASE[randomIndex(pBaseSub)]
		'''
		return (a, y)


	#output the probability of the nucleotide a in sequence x matching the nucleotide b in sequence y
	@staticmethod	
	def emissionProb(a, b, t1, t2):
		f = 0.0
		dN = [a, b]	
		t = [t1, t2]
		for n in settings.BASE:
			fTemp = settings.BASE_FREQUENCY[settings.BASE.index(n)] 
			for i in range(2):
				if n != dN[i]:
					fTemp *= dnaModelJC69.p1(t[i])
				else:
					fTemp *= dnaModelJC69.p0(t[i])
			f += fTemp 	
		for d in dN:
			f = f / settings.BASE_FREQUENCY[settings.BASE.index(d)]
		return f
		
def sampleNucleotide():
        return settings.BASE[randomIndex(settings.BASE_FREQUENCY)]

class xInsertModel(eModel):
	@staticmethod
	def emission(t1, t2, a=None):
		return (sampleNucleotide(), settings.GAP)

	@staticmethod
	def emissionProb(a, b, t1, t2):
		return settings.BASE_FREQUENCY[settings.BASE.index(a)]
		
class yInsertModel(eModel):
	@staticmethod
	def emission(t1, t2, a=None):
		return (settings.GAP, sampleNucleotide())

	@staticmethod
	def emissionProb(a, b, t1, t2):
		return settings.BASE_FREQUENCY[settings.BASE.index(b)]




class State:
	parent = []
	inFreq = []
	children = []
	outFreq = []
	#default dummystate
	stateType = 'D'

	def emission():
		pass
	
	def emissionProb():
		pass

	def __init__(self):
		pass

	def __init__(self, parent, inFreq, children, outFreq, eModel, Type = 'D'):
		self.parent = parent
		self.inFreq = inFreq
		self.children = children
		self.outFreq = outFreq
		self.emission = eModel.emission
		self.emissionProb = eModel.emissionProb
		self.stateType = Type
	
	def nextState(self):
		if 0 == len(self.children):
			print "Reached the end state."
			return False
		else:
			return self.children[randomIndex(self.outFreq)]
#empty for dummystate

def openGapProb(t):
		return (1- math.exp(-settings.INDEL * (t)))

class pHMM:
	gapOpenProb = 1	
	
	xState = None
	yState = None
	mState = None
	wState = None
	time1 = 0
	time2 = 0
	left = "left"
	right = "right"
	rep = 0 #replicate index

	observedDif = 0 #the number of mutation sites
	indelNum = 0 #the number of insertion-deletion events
	matchNum = 0
	#transition Matrix: W X Y M
	transMatrix = []
	stateDistr = None

	def __init__(self,t1, t2, name1="left", name2="right"):
		self.time1 = t1
		self.time2 = t2
		self.left = name1
		self.right = name2
		self.gapOpenProb = openGapProb(t1+t2)
		self.stateInit()
		dnaModelJC69.setEmissionProb(t1+t2)
		#random.seed(1)

	def stateInit(self):
		self.dummyState = State([], [1 - settings.EPSILON, 1- settings.EPSILON, 1 - settings.GAMMA], [None],[1.0], eModel) 
		self.xState = State([], [settings.EPSILON, self.gapOpenProb] ,[None, self.dummyState], [settings.EPSILON ,1-settings.EPSILON], xInsertModel, "X")
		self.xState.children[0] = self.xState
		self.xState.parent.append(self.xState)

		#print "x state out prob:"
		#print self.xState.outFreq

		self.yState = State([], [settings.EPSILON, self.gapOpenProb], [None, self.dummyState], [settings.EPSILON, 1-settings.EPSILON], yInsertModel, "Y")
		self.yState.children[0] = self.yState
		self.yState.parent.append(self.yState)

		#print "y state out prob:"
		#print self.yState.outFreq
	
		self.mState = State([], [settings.GAMMA, 1 - 2 * self.gapOpenProb], [None, self.dummyState], [settings.GAMMA, 1-settings.GAMMA], dnaModelJC69, "M")
		self.mState.children[0] = self.mState
		self.mState.parent.append(self.mState)

		#print "match state out prob:"
		#print self.mState.outFreq

		self.dummyState.parent.append(self.xState)
		self.dummyState.parent.append(self.yState)
		self.dummyState.parent.append(self.mState)
		
		
		self.wState = State([None], [1], [self.xState, self.yState, self.mState],[self.gapOpenProb, self.gapOpenProb, 1-2*self.gapOpenProb], eModel, "W")

		self.dummyState.children[0] = self.wState
		self.wState.parent[0] = self.dummyState
		
		self.xState.parent.append(self.wState)
		self.yState.parent.append(self.wState)
		self.mState.parent.append(self.wState)
		
		self.transMatrix = []
		self.transMatrix.append([0, self.gapOpenProb, self.gapOpenProb, 1 - 2*self.gapOpenProb]) 
		self.transMatrix.append([1-settings.EPSILON, settings.EPSILON, 0, 0])
		self.transMatrix.append([1-settings.EPSILON, settings.EPSILON, 0, 0])
		self.transMatrix.append([1, 0, 0, 0])
		#print "waiting state out prob:"
		#print self.wState.outFreq
	
	def generateSeq(self, rootOrLen):
		children1 = seqDNA("%s_%d:%.4f"%(self.left, self.rep, self.time1))
		children2 = seqDNA("%s_%d:%.4f"%(self.right, self.rep, self.time2))
		self.rep += 1

		seqLen = 0 
		curState = self.wState
		rLen = rootOrLen
		rN = None
		self.observedDif = 0
		self.indelNum = 0
		self.matchNum = 0
		subSiteNum = 0
		stateType = ''
		if isinstance(rootOrLen, seqDNA):
			rLen = rootOrLen.site[rootOrLen.length] #get the last site no. of the sequence
			rN = rootOrLen
			
		siteNum = rLen
		preState = "W"
		while seqLen < rLen:
			site = []
			stateType = curState.stateType
			if "D" != stateType:
				preState = stateType
			if None != rN:
				if settings.GAP != rN.seq[seqLen]:
                        		a = curState.emission(self.time1, self.time2, rN.seq[seqLen])
					curState = curState.nextState()
					site = rN.site[seqLen]
				else:
					a = (settings.GAP, settings.GAP)
			else:
				a = curState.emission(self.time1, self.time2)
				curState = curState.nextState()
                        if None != a:
				seqLen += 1
				if "M" == stateType:
					self.matchNum += 1
					#seqLen += 1
					subSiteNum += 1
					if a[0] != a[1]:
						self.observedDif += 1
				if None != rN:
					if "M" != stateType:
						siteNum += 1  
						children1.insertSeq(a[0], siteNum)
						children2.insertSeq(a[1], siteNum)
						continue
				else:
					children1.insertSeq(a[0], site)
					children2.insertSeq(a[1], site)

			if seqLen < rLen and "W" == preState and ( "X" == curState.stateType or "Y" ==curState.stateType):
				self.indelNum += 1
		self.observedDif = (1.0 * self.observedDif) / subSiteNum
		#print dnaModelJC69
		
		return (children1, children2)	

	@staticmethod	
	def optimalAlignment(seq, score, pointer):
		x = seq[0]
		y = seq[1]
		xLen = len(x.seq)
		yLen = len(y.seq)	
		alignedX = seqDNA(x.name)
		alignedY = seqDNA(y.name)
		pHMM.observedDif = 0
		subSiteNum = 0
		# backtrack start point
		xIndex = xLen
		yIndex = yLen
		optimalState = ""
		optimalScore = 0

		if score[0] > score[1]:
			optimalState = "M"
			optimalScore = score[0]
			if score[0] < score[2]:
				optimalState = "Y"
				optimalScore = score[2]
		else:
			optimalState = "X"
			optimalScore = score[1]
			if score[1] < score[2]:
				optimalState = "Y"
				optimalScore = score[2]
			
		#print "The optimal alignment with log probability %f"%(optimalScore)

		while xIndex > 0 and yIndex > 0:
			nextState = pointer[optimalState][xIndex][yIndex]	
			if "X" == optimalState:
				xIndex -= 1
				alignedX.insertSeq(x.seq[xIndex], x.site[xIndex])
				alignedY.insertSeq('-', 0)
			elif "Y" == optimalState:
				yIndex -= 1
				alignedX.insertSeq('-',0)
				alignedY.insertSeq(y.seq[yIndex], y.site[yIndex])
			elif "M" == optimalState:
				xIndex -= 1
				yIndex -=1
				alignedX.insertSeq(x.seq[xIndex], x.site[xIndex])
				alignedY.insertSeq(y.seq[yIndex],y.site[yIndex])
				subSiteNum += 1
				if x.seq[xIndex] != y.seq[yIndex]:
					pHMM.observedDif += 1
			optimalState = nextState

		if xIndex > 0:
			alignedX.insertSeq(x.seq[:xIndex][::-1],x.site[:xIndex][::-1])
			alignedY.insertSeq('-'*xIndex, [0 for i in range(xIndex)])
		elif yIndex > 0:
			alignedY.insertSeq(y.seq[:yIndex][::-1], y.site[:yIndex][::-1])
			alignedX.insertSeq('-'*yIndex, [0 for i in range(yIndex)])

		pHMM.observedDif = (1.0 * pHMM.observedDif)/subSiteNum
		alignedX.seq = alignedX.seq[::-1]
		alignedY.seq = alignedY.seq[::-1]
		alignedX.site = alignedX.site[::-1]
		alignedY.site = alignedY.site[::-1]

		return (alignedX, alignedY)
					

			
		
	
	def alignSeq(self, seq):
		eProb = lambda state,x,y: state.emissionProb(x,y,self.time1, self.time2)
		seqX = seq[0].seq
		seqY = seq[1].seq

		if len(seqX) < 1 or len(seqY) < 1:
			print "At least one sequence is empty."
			sys.exit(0)

		xRange = len(seqX) + 1
		yRange = len(seqY) + 1

		def log(x):
			if x <= 0:
				return float("-inf")
			else:
				return math.log(x)	

		'''
			Initial the matrix for score and pointer used to traceback the path.
		'''
		mInf = float('-Inf')
		wMatrix = [mInf for j in range(yRange)]
		scoreMatrix = {"M": copy.deepcopy(wMatrix), "X": copy.deepcopy(wMatrix), "Y": copy.deepcopy(wMatrix)}

		pMatrix = [[ 0 for j in range(yRange)] for i in range(xRange)]
		pointerMatrix = {"M": copy.deepcopy(pMatrix), "X": copy.deepcopy(pMatrix), "Y": copy.deepcopy(pMatrix)}

		# The initial value 
		scoreMatrix["M"][0] = 0
		pointerMatrix["M"][0][0] = "start"


		stateArr = [self.mState, self.xState, self.yState]

		def getValue(scoreMatrix,stateType,j):
			if j < 0:
				return float("-Inf")
			else:
				return scoreMatrix[stateType][j]
		stateNum = len(stateArr)
		a = [mInf for i in range(stateNum)]
		s = [0 for i in range(stateNum)]
		start = time.clock()
		for i in range(xRange):
			for j in range(yRange):
				if i != 0 or j != 0:
					maxScore = mInf
					#update the score and pointer matrix
					for state in stateArr:
						yIndex = j - 1

						if "X" == state.stateType:
							yIndex = j
						maxPreScore = getValue(scoreMatrix, state.stateType, yIndex) + log(state.inFreq[0])
						tempPointer = state.stateType
						for stateIndex in range(len(self.dummyState.parent)):
							stateType = self.dummyState.parent[stateIndex].stateType
							tempProb = getValue(scoreMatrix, stateType, yIndex) + log(self.dummyState.inFreq[stateIndex]) + log(state.inFreq[1])
							if tempProb > maxPreScore:
								tempPointer = stateType
								maxPreScore = tempProb

						pointerMatrix[state.stateType][i][j] = tempPointer
						
						if (i > 0 and j > 0) or (i > 0 and "X" == state.stateType) or (j > 0 and "Y" == state.stateType):		
							s[stateArr.index(state)] = log(eProb(state, seqX[i-1],seqY[j-1])) + maxPreScore
						else:
							s[stateArr.index(state)] = mInf + maxPreScore 	

						if "M" == state.stateType and j - 1 > 0:
							for stateIndex in range(stateNum):
								scoreMatrix[stateArr[stateIndex].stateType][yIndex] = a[stateIndex]
					if xRange - 1 != i or yRange - 1 != j:
						a = copy.deepcopy(s)
		return self.optimalAlignment(seq, s, pointerMatrix)

	def setStateDistribution(self):
		np.savetxt('./data/input', np.array(self.transMatrix))
		os.system("./stateDistribution.o>./data/result;")
		rf = open('./data/result', 'r')
		s = rf.read()
		rf.close()
		self.stateDistr = [float(i) for i in s.strip('\n').split()]

	def getExpectedInDel(self):
		#print "transition matrix w x y m(t:%.4f):"%(settings.TIME)
		#print self.transMatrix
		self.setStateDistribution()
		#print "state distribution:"
		#print self.stateDistr

		wStateNum = settings.LENGTH / (sum(self.stateDistr[1:])) * self.stateDistr[0]
		return (wStateNum * sum(self.transMatrix[0][1:3]))

	def getExpectedMatch(self):
		if None == self.stateDistr:
			self.setStateDistribution()
		return (settings.LENGTH * (self.stateDistr[3] / sum(self.stateDistr[1:])))

