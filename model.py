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
	def emissionProb(a,b,t1, t2):#log probability
		return float("-Inf")

		
class dnaModelJC69(eModel):

	modelName = "JC69"	
	matchProb = 0
	mismatchProb = 0
	# remain probability	
	@staticmethod
	def p0(t):
		return (1.0/4 + 3.0/4 * math.exp(-4 * settings.LAMBDA * (t/3.0))) 

	# substitution probability
	@staticmethod
	def p1(t):
		return (1.0/4 - 1.0/4 * math.exp(-4 * settings.LAMBDA * (t/3.0)))

	eMatrix = None	
	
	@staticmethod
	def setEmissionProb(t):
		dnaModelJC69.eMatrix = [[dnaModelJC69.p1(t) for i in range(4)] for j in range(4)]
		for i in range(4):
			dnaModelJC69.eMatrix[i][i] = dnaModelJC69.p0(t)
		'''
		dnaModelJC69.matchProb = math.log(0.25) + math.log(3*math.pow(dnaModelJC69.p1(t/2.0), 2) + math.pow(dnaModelJC69.p0(t/2.0), 2))# - math.log(math.pow(0.25, 2))
		dnaModelJC69.mismatchProb = math.log(0.5) + math.log(dnaModelJC69.p1(t/2.0)) + math.log(dnaModelJC69.p1(t/2.0)+dnaModelJC69.p0(t/2.0))# - math.log(math.pow(0.25, 2))
		'''

	
	@staticmethod
	def emission(t1, t2, a=None):
		if settings.GAP == a:
			return (settings.GAP, settings.GAP)
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
		if a == b:
			return math.log(dnaModelJC69.p0(t1+t2)) + math.log(settings.BASE_FREQUENCY[settings.BASE.index(a)])
		else:
			return math.log(dnaModelJC69.p1(t1+t2)) + math.log(settings.BASE_FREQUENCY[settings.BASE.index(a)])
		
def sampleNucleotide():
        return settings.BASE[randomIndex(settings.BASE_FREQUENCY)]

class xInsertModel(eModel):
	@staticmethod
	def emission(t1, t2, a=None):
		return (sampleNucleotide(), settings.GAP)

	@staticmethod
	def emissionProb(a, b, t1, t2):
		return math.log(settings.BASE_FREQUENCY[settings.BASE.index(a)])
		
class yInsertModel(eModel):
	@staticmethod
	def emission(t1, t2, a=None):
		return (settings.GAP, sampleNucleotide())

	@staticmethod
	def emissionProb(a, b, t1, t2):
		return math.log(settings.BASE_FREQUENCY[settings.BASE.index(b)])




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
	#extension:{indelRate:{t:openGapProb}}
	rateToProb={'0.4':{'0.1':{'0.1':0.05, '0.2':0.01, '0.3':0.015,'0.4':0.02,'0.5':0.025,'0.6':0.03},
			 '0.2':{'0.1':0.01, '0.2':0.02, '0.3':0.03,'0.4':0.04,'0.5':0.05,'0.6':0.06},
			 '0.3':{'0.4':0.06},
			 '0.4':{'0.4':0.08},
			 '0.5':{'0.4':0.1},
			 '0.6':{'0.4':0.12}}}

	def __init__(self,t1, t2, name1="left", name2="right"):
		self.time1 = t1
		self.time2 = t2
		self.left = name1
		self.right = name2
		self.gapOpenProb = self.openGapProb(t1+t2)
		self.stateInit()
		dnaModelJC69.setEmissionProb(t1+t2)
		#random.seed(1)
	
	def openGapProb(self, t): 
		return self.rateToProb['%.1f'%settings.EPSILON]['%.1f'%settings.INDEL]['%.1f'%t]
	def stateInit(self):
		self.xState = State([], [settings.EPSILON, self.gapOpenProb, self.gapOpenProb] ,[], [settings.EPSILON ,1-settings.EPSILON], xInsertModel, "X")
		self.xState.children.append(self.xState)
		self.xState.parent.append(self.xState)


		self.yState = State([], [settings.EPSILON, self.gapOpenProb, self.gapOpenProb], [], [settings.EPSILON, 1-settings.EPSILON], yInsertModel, "Y")
		self.yState.children.append(self.yState)
		self.yState.parent.append(self.yState)

	
		self.mState = State([None, self.xState, self.yState], [1 - 2 * self.gapOpenProb, 1 - settings.EPSILON, 1 - settings.EPSILON, 1 - 2 * self.gapOpenProb], [None, self.xState, self.yState], [1 - 2 * self.gapOpenProb, self.gapOpenProb, self.gapOpenProb], dnaModelJC69, "M")
		self.mState.children[0] = self.mState
		self.mState.parent[0] = self.mState
	
		self.wState = State([], [], [self.xState, self.yState, self.mState], [self.gapOpenProb, self.gapOpenProb, 1 - 2 * self.gapOpenProb], eModel, "W")
		self.mState.parent.append(self.wState)
		self.xState.children.append(self.mState)
		self.yState.children.append(self.mState)

		self.xState.parent.append(self.mState)
		self.xState.parent.append(self.wState)
		self.yState.parent.append(self.mState)
		self.yState.parent.append(self.wState)
		
		self.transMatrix = []
		self.transMatrix.append([0, self.gapOpenProb, self.gapOpenProb, 1 - 2 * self.gapOpenProb])
		self.transMatrix.append([0, settings.EPSILON, 0, 1 - settings.EPSILON])
		self.transMatrix.append([0, 0, settings.EPSILON, 1 - settings.EPSILON])
		self.transMatrix.append([0, self.gapOpenProb, self.gapOpenProb, 1 - 2 * self.gapOpenProb])

	
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
		a = None
		while seqLen < rLen:
			site = []
			stateType = curState.stateType

			if None != rN:
				a = curState.emission(self.time1, self.time2, rN.seq[seqLen])
			else:
				a = curState.emission(self.time1, self.time2)

			nextState = curState.nextState()
			if None != a:
				seqLen += 1
				if None != rN:
					if settings.GAP != rN.seq[seqLen] and "M" != stateType:
						siteNum += 1
						children1.insertSeq(a[0], siteNum)
						children2.insertSeq(a[1], siteNum)
				else:
					children1.insertSeq(a[0], site) 
					children2.insertSeq(a[1], site) 

				if "M" == stateType:
					self.matchNum += 1
					if a[0] != a[1]:
						self.observedDif += 1

			if seqLen < rLen and nextState != curState and ("X" == nextState.stateType or "Y" == nextState.stateType):
				self.indelNum += 1
			curState = nextState

		self.observedDif = (1.0 * self.observedDif) / self.matchNum
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
			optimalState = "X"
			optimalScore = score[0]
			if score[0] < score[2]:
				optimalState = "M"
				optimalScore = score[2]
		else:
			optimalState = "Y"
			optimalScore = score[1]
			if score[1] < score[2]:
				optimalState = "M"
				optimalScore = score[2]
		#print "The optimal alignment with log probability %f"%(optimalScore)

		while xIndex > 0 and yIndex > 0:
			#print xIndex,yIndex,optimalState
			nextState = pointer[optimalState][xIndex][yIndex]	
			#print nextState
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
				yIndex -= 1
				alignedX.insertSeq(x.seq[xIndex], x.site[xIndex])
				alignedY.insertSeq(y.seq[yIndex], y.site[yIndex])
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
		if 0 != subSiteNum:
			pHMM.observedDif = (1.0 * pHMM.observedDif)/subSiteNum

		alignedX.seq = alignedX.seq[::-1]
		alignedY.seq = alignedY.seq[::-1]
		alignedX.site = alignedX.site[::-1]
		alignedY.site = alignedY.site[::-1]

		return ([alignedX, alignedY], optimalScore)
					

			
		
	
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
		wMatrix = [[mInf for j in range(yRange)] for i in range(xRange)]
		scoreMatrix = {"X": copy.deepcopy(wMatrix), "Y": copy.deepcopy(wMatrix), "M": copy.deepcopy(wMatrix)}

		pMatrix = [[ 0 for j in range(yRange)] for i in range(xRange)]
		pointerMatrix = {"X": copy.deepcopy(pMatrix), "Y": copy.deepcopy(pMatrix), "M": copy.deepcopy(pMatrix)}

		# The initial value 
		scoreMatrix["M"][0][0] = 0
		pointerMatrix["M"][0][0] = "start"


		stateArr = [self.xState, self.yState, self.mState]


		for i in range(1, xRange):
			for j in range(1, yRange):
				for state in stateArr:
					xIndex = i - 1
					yIndex = j - 1
					if "X" == state.stateType:
						yIndex += 1
					elif "Y" == state.stateType:
						xIndex += 1
					maxScore = scoreMatrix[state.stateType][xIndex][yIndex] + log(state.inFreq[0])
					tempPointer = state.stateType
					for stateIndex in range(1, len(state.parent)-1):
						stateType = state.parent[stateIndex].stateType
						tempScore = scoreMatrix[stateType][xIndex][yIndex] + log(state.inFreq[stateIndex])
						if tempScore > maxScore:
							maxScore = tempScore
							tempPointer = stateType
				
					scoreMatrix[state.stateType][i][j] = eProb(state, seqX[i-1], seqY[j-1]) + maxScore 
					pointerMatrix[state.stateType][i][j] = tempPointer
		'''
		for p in pointerMatrix:
			print p
			for i in range(xRange):
				print pointerMatrix[p][i]
		'''
		s = []
		for state in stateArr:
			s.append(scoreMatrix[state.stateType][xRange-1][yRange-1])
		return self.optimalAlignment(seq, s, pointerMatrix)
	
	def getProbability(self, alignment):
		eProb = lambda state,a,b: state.emissionProb(a, b, self.time1, self.time2)
		x = alignment[0].seq
		y = alignment[1].seq
		score = 0
		preState = 0
		stateArr = [self.wState, self.xState, self.yState, self.mState]
		
		for index in range(alignment[0].length):
			curState = 3
			if settings.GAP == x[index]:#yState
				curState = 2
			elif settings.GAP == y[index]:#xState 
				curState = 1
			score += eProb(stateArr[curState], x[index], y[index])
			if 0 == index:
				score += math.log(stateArr[preState].outFreq[curState - 1])	
			else:
				if preState == curState:
					score += math.log(stateArr[preState].outFreq[0])
				else:
					if 3 == preState:
						score += math.log(stateArr[preState].outFreq[curState])
					else:
						score += math.log(stateArr[preState].outFreq[1])
			preState = curState

		return score

	def setStateDistribution(self):
		np.savetxt('input', np.array(self.transMatrix))
		os.system("./stateDistribution.o>./data/result;")
		rf = open('./data/result', 'r')
		s = rf.read()
		rf.close()
		self.stateDistr = [float(i) for i in s.strip('\n').split()]

	def getExpectedInDel(self):
		self.setStateDistribution()
		mStateNum = settings.LENGTH / sum(self.stateDistr[1:]) * self.stateDistr[3]
		wStateNum = settings.LENGTH / sum(self.stateDistr[1:]) * self.stateDistr[0]
		return wStateNum * sum(self.transMatrix[0][1:3]) + mStateNum * sum(self.transMatrix[3][1:3])

	def getExpectedMatch(self):
		if None == self.stateDistr:
			self.setStateDistribution()
		return  settings.LENGTH * (self.stateDistr[3] / sum(self.stateDistr[1:]))

class PAGAN(pHMM):
	
	def openGapProb(self, t): 
		return (1- math.exp(-settings.INDEL * (t))) 

	def stateInit(self):
		self.dummyState = State([], [1 - settings.EPSILON, 1- settings.EPSILON, 1 - settings.GAMMA], [None],[1.0], eModel) 
		self.xState = State([], [settings.EPSILON, self.gapOpenProb] ,[None, self.dummyState], [settings.EPSILON ,1-settings.EPSILON], xInsertModel, "X")
		self.xState.children[0] = self.xState
		self.xState.parent.append(self.xState)


		self.yState = State([], [settings.EPSILON, self.gapOpenProb], [None, self.dummyState], [settings.EPSILON, 1-settings.EPSILON], yInsertModel, "Y")
		self.yState.children[0] = self.yState
		self.yState.parent.append(self.yState)

	
		self.mState = State([], [settings.GAMMA, 1 - 2 * self.gapOpenProb], [None, self.dummyState], [settings.GAMMA, 1-settings.GAMMA], dnaModelJC69, "M")
		self.mState.children[0] = self.mState
		self.mState.parent.append(self.mState)


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
		self.transMatrix.append([1-settings.EPSILON, 0, settings.EPSILON, 0])
		self.transMatrix.append([1, 0, 0, 0])

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
			optimalState = "X"
			optimalScore = score[0]
			if score[0] < score[2]:
				optimalState = "M"
				optimalScore = score[2]
		else:
			optimalState = "Y"
			optimalScore = score[1]
			if score[1] < score[2]:
				optimalState = "M"
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
		if 0 != subSiteNum:
			pHMM.observedDif = (1.0 * pHMM.observedDif)/subSiteNum
		alignedX.seq = alignedX.seq[::-1]
		alignedY.seq = alignedY.seq[::-1]
		alignedX.site = alignedX.site[::-1]
		alignedY.site = alignedY.site[::-1]

		return ([alignedX, alignedY], optimalScore)
		
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
		wMatrix = [[mInf for j in range(yRange)] for i in range(xRange)]
		scoreMatrix = {"X": copy.deepcopy(wMatrix), "Y": copy.deepcopy(wMatrix), "M": copy.deepcopy(wMatrix)}

		pMatrix = [[ 0 for j in range(yRange)] for i in range(xRange)]
		pointerMatrix = {"X": copy.deepcopy(pMatrix), "Y": copy.deepcopy(pMatrix), "M": copy.deepcopy(pMatrix)}

		# The initial value
		scoreMatrix["M"][0][0] = 0
		pointerMatrix["M"][0][0] = "start"


		stateArr = [self.xState, self.yState, self.mState]
		
		for i in range(1, xRange):
			for j in range(1, yRange):
				for state in stateArr:
					xIndex = i - 1
					yIndex = j - 1
					if "X" == state.stateType:
						yIndex += 1
					elif "Y" == state.stateType:
						xIndex += 1
					maxScore = scoreMatrix[state.stateType][xIndex][yIndex] + log(state.inFreq[0])
					tempPointer = state.stateType
					for stateIndex in range(len(self.dummyState.parent)):
						stateType = self.dummyState.parent[stateIndex].stateType
						tempScore = scoreMatrix[stateType][xIndex][yIndex] + log(self.dummyState.inFreq[stateIndex]) + log(state.inFreq[1])
						if tempScore > maxScore:
							maxScore = tempScore
							tempPointer = stateType
					scoreMatrix[state.stateType][i][j] = eProb(state, seqX[i-1], seqY[j-1]) + maxScore
					pointerMatrix[state.stateType][i][j] = tempPointer
		s = []
		for state in stateArr:
			s.append(scoreMatrix[state.stateType][xRange - 1][yRange - 1])
			
		return self.optimalAlignment(seq, s, pointerMatrix)

	def getProbability(self, alignment):
		eProb = lambda state,a,b: state.emissionProb(a, b, self.time1, self.time2)
		x = alignment[0].seq
		y = alignment[1].seq
		score = 0
		preState = 0
		stateArr = [self.xState, self.yState, self.mState]
		
		for index in range(alignment[0].length):
			curState = 2
			if settings.GAP == x[index]:
				curState = 1
			elif settings.GAP == y[index]: 
				curState = 0
			score += eProb(stateArr[curState], x[index], y[index])
			if 0 == index:
				score += math.log(self.wState.outFreq[curState])
			else:
				if preState != curState:
					score += math.log(self.dummyState.inFreq[preState])
					score += math.log(self.wState.outFreq[curState])
				elif 2 == curState:
					score += math.log(self.wState.outFreq[curState])
				else:
					score += math.log(stateArr[curState].inFreq[0])
			preState = curState
		return score

	def getExpectedInDel(self):
		self.setStateDistribution()
		wStateNum = settings.LENGTH / (sum(self.stateDistr[1:])) * self.stateDistr[0]
		return (wStateNum * sum(self.transMatrix[0][1:3]))

