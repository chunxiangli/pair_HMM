#!/usr/bin/python
import os
import sys
import settings
import re
import random
import copy
import math
import time
import multiprocessing
from input import *
from model import *
from utils.myprint import *
from utils.tools import *
from utils.myplot import *
from multiprocessing import Pool
from Bio import SeqIO
import multiprocessing as mp
import numpy as np
from matplotlib import pyplot
#from scipy import stats


cpuNum  = multiprocessing.cpu_count()
sequencesBlock = []
def printHelp(argv):
	print '''
	Usage: %s [args]
	Example: %s len=20 
	
	You can either edit settings.py to handle algorithms' setting or use command line argument (key=value). Command line argument will override the setting defined in settings.py when run.

	General command line arguments are:
		len : The length of randomly generated ancestral sequence
		  t : The evolution time difference between the two generated sequences
	       step : Increment 
	       iter : 
		  l : Lambda, the parameter for Jc69 model. The probability for one nucleotide substituted by another nucleotides is 1/4 -(1/4)exp(-4lt), remain is 1/4 + (3/4)exp(-4lt) 
		  g : Gamma, the match extension probability
		  e : gap extension probability
		  a : indel rate. The gap open probability is 1 - exp(-at)
		 in : input file or the output file for tree simulation result
		out : output the alignment result
	       tree : tree architecure discription. i.e. (t1:0.1,t2:0.1)
	'''%(argv[0],argv[0]) 
	

def validateArgv(argv):
	for x in range(1, len(argv)):
		arg = argv[x].split("=", 2)
		
		if arg[0] == 'help':
			printHelp(argv)
			sys.exit(0)

		if len(arg) != 2:
			print "Argument \'%s\' in wrong format. Use \'key=value\'." % (arg[0])
			print "For more instructions try \'%s help\'" % (argv[0])
			sys.exit(0)	

def poolMap(task, param, cpuNum, chuckSize=1):
	pool = Pool(cpuNum)
	result = pool.map(task, param, chuckSize)
	pool.close()
	pool.join()
	return result

def fidelity(real, path):
	rNum = 0
	assert 2 == len(real), "Invalide real alignment."
	assert 2 == len(path), "Invalide alignment."
	'''
	print "original"
	printPairAlign([real[0].seq,real[1].seq])
	print "alignment"
	printPairAlign([path[0].seq,path[1].seq])
	'''
	for index in range(real[0].length):
		a = real[0].seq[index]
		b = real[1].seq[index]
		if not (settings.GAP == a and settings.GAP == b):
			site = real[0].site[index]
			if settings.GAP == a:
				pIndex = path[1].site.index(site)
				if settings.GAP == path[0].seq[pIndex]:
					rNum += 1
			elif settings.GAP == b:
				pIndex = path[0].site.index(site)
				if settings.GAP == path[1].seq[pIndex]:
					rNum += 1
			else:
				pIndex = path[0].site.index(site)
				if path[1].site[pIndex] == site:
					rNum += 1
	return ((1.0*rNum)/real[0].length)

def observedDifference(seqPair):
	assert 2 == len(seqPair), "Please give a pair alignment sequences."
	assert isinstance(seqPair[0], seqDNA), "The format of the first sequence wrong."
	assert isinstance(seqPair[1], seqDNA), "The format of the second sequence wrong."
	observedP = 0
	subSiteNum = 0
	for index in range(seqPair[0].length):
		if settings.GAP != seqPair[0].seq[index] and settings.GAP != seqPair[1].seq[index]:
			subSiteNum += 1
			if seqPair[0].seq[index] != seqPair[1].seq[index]:
				observedP += 1 
	observedP = (1.0 * observedP) / subSiteNum
	if observedP >= 3.0/4:
		print "Can't estimate the distance since the observed difference is too large."
		return (-1)
	return (observedP) 

	
def pairSimulation(outFile, plot=False):
	if None == settings.TREE:
		(left, t1) = ('A', (settings.TIME/2))
		(right, t2) = ('B', (settings.TIME/2))
	else:
		(left, t1) = findTaxaAndLength(settings.TREE)	
		(right, t2) = findTaxaAndLength(settings.TREE[len(left)+len(t1)+2:])
		t1 = float(t1)
		t2 = float(t2)
		settings.TIME = t1 + t2
	pair_hmm = None
	if "Durbin" == settings.MODEL:
		pair_hmm = pHMM(t1, t2, left, right)
	elif "aModel" == settings.MODEL:
		pair_hmm = PAGAN(t1, t2, left, right)
	#np.savetxt('%s/substitutionMatrix'%(os.path.dirname(settings.IN_FILE)), np.array(dnaModelJC69.eMatrix), fmt="%.5f")	
	p = []
	indelArr = []
	matchArr = []
	for i in range(settings.REPLICATE):
		children  = pair_hmm.generateSeq(settings.LENGTH)
		printFastaFile(children, outFile)
                if pair_hmm.observedDif >= 0.75:
			print "To large difference"
			p.append(-1)
		else:
			p.append(pair_hmm.observedDif)
                indelArr.append(pair_hmm.indelNum)
                matchArr.append(pair_hmm.matchNum)
	if plot and settings.REPLICATE > 2:
		eInDel = pair_hmm.getExpectedInDel()
		eMatch = pair_hmm.getExpectedMatch()
		fName = "estimateDistance_%s.png"%(re.sub(r'\.fas', '', os.path.basename(settings.IN_FILE)))
		xlabel = "Replicate Sequence"
		title = "Simulation(a=%.3f,e=%.3f,l=%.2f,g=%.2f)"%(settings.INDEL, settings.EPSILON, settings.LAMBDA, settings.GAMMA)
        	plotSimulationResult(np.arange(settings.REPLICATE), t1+t2, d2t(p2d(np.array(p))), None, eInDel, np.array(indelArr), eMatch, np.array(matchArr), xlabel, fName, title)

def readPairSequences(inFile, num=1):
	pairArr = []
	sequences = []
	def name(x):
		if settings.REAL_TIME:
			return "%s:%.2f"%(x, settings.REAL_TIME)
		return x
			
	for s in SeqIO.parse(inFile, settings.IN_FORMAT):
		if len(pairArr) < num:
			sequences.append(seqDNA(name(s.id), s.seq.data))
			if 2 == len(sequences):
				pairArr.append(sequences[:])
				sequences = []
	return pairArr
			
def pairAlignmentWithSpecificParameters(param):
	t1 = settings.TIME/2
	t2 = settings.TIME/2
	sequences =  param[0]
	realTime1 = float(re.search(':(\d+\.?\d+)\|?', sequences[0].name).group(1).strip())
	realTime2 = float(re.search(':(\d+\.?\d+)\|?', sequences[1].name).group(1).strip()) 
	if not param[1]:
		t1 = realTime1
		t2 = realTime2
	pair_hmm = None
	if "Durbin" == settings.MODEL:
		pair_hmm = pHMM(t1, t2)
	elif "aModel" == settings.MODEL:
		pair_hmm = PAGAN(t1, t2)
	alignment = copy.deepcopy(sequences)
	for c in alignment: c.removeGap()
	(alignment, score) = pair_hmm.alignSeq(alignment)	
	if param[1]:
		pHmm = None
		if settings.ACHECK:
			settings.INDEL = float(re.search('_a(\d+\.?\d+)_', settings.IN_FILE).group(1))
		if "Durbin" == settings.MODEL:
			pHmm = pHMM(realTime1, realTime2) 
		elif "aModel" == settings.MODEL:
			pHmm = PAGAN(realTime1, realTime2)
		score = pHmm.getProbability(alignment)
	accuracy = fidelity(sequences, alignment)
	#np.savetxt("%s/substitutionMatrix_%.2f"%(os.path.dirname(settings.IN_FILE), settings.TIME), np.array(dnaModelJC69.eMatrix), fmt="%.5f")
	return [accuracy, alignment, score]

def chunkSize(size):
	cSize = size / cpuNum
	if cSize * cpuNum < size:
		cSize += 1
	return cSize

def realignment(inFile, outFile, reset=False):
	sequences= readPairSequences(inFile, settings.REPLICATE)
	pArr = [[seq, reset] for seq in sequences]
	results = poolMap(pairAlignmentWithSpecificParameters, pArr, chunkSize(settings.REPLICATE))
	results = np.array(results)
	printFastaFile(results[:,1], outFile)
	return (results[:,0], results[:,2])
	
	
def alignmentWithSpecificDistance(t):
        f1 = open(settings.IN_FILE, 'r')
        settings.TIME = t
        f2 = open("%s/alignment_over_%s"%(os.path.dirname(settings.OUT_FILE), re.sub(r'\.fas$', '_with_a%.2f_and_t%.2f.fas'%(settings.INDEL, settings.TIME), os.path.basename(settings.IN_FILE))),'w')
        (accArr, scoreArr) = realignment(f1, f2, True)
	#print accArr
	#print scoreArr
        f2.close()
        f1.close()
        return (accArr.mean(), scoreArr.mean())

def simulationAndAlignmentWithVariantDistance(plot=False):
	P = []
	t = [ 0.05+i*0.01 for i in range(36)]
	accuracy = []
	realTime = 0
	try:
		replicateSimulation(settings.IN_FILE)
		indel = settings.INDEL
		realTime = settings.TIME
		start = time.time()
		for i in t:
			accuracy.append(alignmentWithSpecificDistance(i)[0])
		print "time using:%.2gs"%(time.time() - start)
		
		if plot:
			tmpS = "(a=%.3f,e=%.3f,l=%.2f,g=%.2f, len=%d, rep=%d)"%(settings.INDEL, settings.EPSILON, settings.LAMBDA, settings.GAMMA, settings.LENGTH, settings.REPLICATE)
			fName = "%s_alignmentWithVariantDistance_over_%s.png"%(settings.MODEL, re.sub(r'\.fas', '', os.path.basename(settings.IN_FILE)))
			title = "Alignment%s"%(tmpS)
			plotAlignmentResult(t, [accuracy], realTime ,fName, title)

	except IOError as e:
        	print 'I/0 error{0}:{1}'.format(e.errno, e.strerror)
		sys.exit(0)
	return accuracy

def pairSimulationAndAlignment(plot=False):
	try:
		replicateSimulation(settings.IN_FILE)
		f1 = open(settings.IN_FILE, 'r')
        	f2 = open("%s/alignment_over_%s"%(os.path.dirname(settings.IN_FILE), os.path.basename(settings.IN_FILE)), 'w')
        	accArr = realignment(f1, f2)
        	f2.close()
        	f1.close()
		if plot and settings.REPLICATE > 3:
                        fName = "pairAlignment_over_%s.png"%(re.sub(r'\.fas', '', os.path.basename(settings.IN_FILE)))
                        title = "Alignment(a=%.2f,e=%.2f,l=%.2f,g=%.2f, len=%d, rep=%d)"%(settings.INDEL, settings.EPSILON, settings.LAMBDA, settings.GAMMA, settings.LENGTH, settings.REPLICATE)
                        plotAlignmentResult(np.arange(settings.REPLICATE), [accArr], settings.TIME, fName, title)
	except IOError as e:
		print 'I/0 error{0}:{1}'.format(e.errno, e.strerror)
                sys.exit(0)

	return (accArr)

def replicateSimulation(outFile, plot=False):
	of = open(outFile, 'w')
	pairSimulation(of, plot)
	of.close()	

def checkAlignmentParameters():
	indel = settings.INDEL
	replicateSimulation(settings.IN_FILE, False)
	realTime = settings.TIME
	
	indelArr = [ 0.02 + i * 0.01 for i in range(9)]
	tArr = [0.05 + i * 0.01 for i in range(36)]
	accuracy = []
	start = time.time()
	for i in indelArr:
		settings.INDEL = i
		acc = []
		for j in tArr:	
			acc = alignmentWithSpecificDistance(j)[0]
		accuracy.append(acc[:])
	print "time using:%.3gs"%(time.time()-start)
	tmpS = "(t%.3f,a%.3f,e%.3f,l%.3f,g%.3f,len%d,rep%d)"%(realTime, indel, settings.EPSILON, settings.LAMBDA, settings.GAMMA, settings.LENGTH, settings.REPLICATE)
	fName = "%s_alignmentOverVariantParameters_over_%s.png"%(settings.MODEL, e.sub(r'\.fas', '', os.path.basename(settings.IN_FILE)))
	title = "Alignment%s"%(tmpS)
	plotAlignmentResult(tArr, accuracy, realTime ,fName, title, indelArr)

def runIndelible(start, step, times, subRate, length):
	observedDif = []
	os.system("./runIndelible %f %f %d %f %d>>log.txt"%(start, step, times, subRate, length))
	for i in range(times):
		iF = "%s/output%.3f.fas"%(settings.IndelibleDataDir, (start + i*step)/2)
		inFile = open(iF, 'r')
		s = inFile.readline()
		sequences = []
		for j in range(2):
			while '' != s and '>' != s[0]: s = inFile.readline()
                	if '' == s:
                        	print "There isn\'t sequence in the fasta file %s."%(iF)
                        	inFile.close()
                        	sys.exit(0)
                	name = s.strip('\n')
                	seq = ''
                	s = inFile.readline()
                	while '' != s and '>' != s[0]:
                        	seq += s
                        	s = inFile.readline()
                	if '' == s and '' == seq:
                        	print "There is only one sequence in the fasta file %s."%(iF)
                        	inFile.close()
                        	sys.exit()
                	sequences.append(seqDNA(name, re.sub(r'\n', '', seq).strip(' ')))
		inFile.close()
		observedDif.append(observedDifference(sequences))

		#printPairAlign([sequences[0].seq, sequences[1].seq])	
	return observedDif

def runPhylosim(start, step, times, subRate, length):
	iF = "phylosim.fas"
	observedDif = []
	os.system("Rscript ./runPhylosim.r %f %f %d %d %s>>log.txt"%(start, step, times, length, iF))
	iF = "./data/%s"%(iF)
        inFile = open(iF, 'r')
        s = inFile.readline()
	for i in range(times):
                sequences = []
                for j in range(2):
                        while '' != s and '>' != s[0]: s = inFile.readline()
                        if '' == s:
                                print "There isn\'t sequence in the fasta file %s."%(iF)
                                inFile.close()
                                sys.exit(0)
                        name = s.strip('\n')
                        seq = ''
                        s = inFile.readline()
                        while '' != s and '>' != s[0]:
                                seq += s
                                s = inFile.readline()
                        if '' == s and '' == seq:
                                print "There is only one sequence in the fasta file %s."%(iF)
                                inFile.close()
                                sys.exit()
                        sequences.append(seqDNA(name, re.sub(r'\n', '', seq).strip(' ')))
                observedDif.append(observedDifference(sequences))
	inFile.close()
	return observedDif

def p2d(x):
	return (-(3.0/4) * np.log(1 - 4/3.0 * x))

def d2t(x):
	return (x/(settings.LAMBDA))

def checkAccuracySimulation():
	t = []
	d = []	
	indel = []
	eInDel = []
	match = []
	eMatch = []
	accuracy = []
	meanEstimateDistance = 0
	meanAccuracy = 0
	localIter = settings.REPLICATE 
	try:
		start = settings.TIME
		#f1 = open(settings.IN_FILE, 'w')
		for i in range(settings.ITERATE):
			t.append(settings.TIME)
			sumD = 0
			wIter = 0
			indelArr = []
			matchArr = []

			pair_hmm = pHMM(settings.TIME/2, settings.TIME/2, "A", "B")
			
			chunkSize = localIter / (5 * 16)	
			if chunkSize * 5 * 16 < localIter:
				chunkSize += 1
			s = 0
			for k in range(chunkSize):
				f1 = open(re.sub(r'\.fas', '_t%.2f_%d.fas'%(settings.TIME, k), settings.IN_FILE), 'w')
				for j in range(80):
					if s < localIter:
        					children  = pair_hmm.generateSeq(settings.LENGTH)
        					printFastaFile(children, f1)
						p = pair_hmm.observedDif
						if p >= 0.75:
							wIter += 1
						else:
							sumD += p2d(p)	
						indelArr.append(pair_hmm.indelNum)
						matchArr.append(pair_hmm.matchNum)
					s += 1
				f1.close()

			eInDel.append(pair_hmm.getExpectedInDel())
			indel.append(np.array(indelArr).mean())
			match.append(np.array(matchArr).mean())
			eMatch.append(pair_hmm.getExpectedMatch())
			d.append(sumD/(localIter - wIter))
			settings.TIME += settings.STEP		
		indelDif = runIndelible(start, settings.STEP, settings.ITERATE, settings.INDEL, settings.LENGTH)
		#phylosimDif = runPhylosim(start, settings.STEP, settings.ITERATE, settings.INDEL, settings.LENGTH)
      		
		if 1 != settings.ITERATE:
			t = np.array(t)
			oDReal = np.array(d)
			oTReal = d2t(oDReal)
			oInDel = np.array(indel)
			oIndelibleDif = d2t(p2d(np.array(indelDif)))
			#oPhylosimDif = d2t(p2d(np.array(phylosimDif)))
			#s = stats.linregress(t, oTReal)
			#lTReal = s[0]*t+s[1]
			fName = "estimate_distance_over_%s.png"%(re.sub(r'\.fas', '', os.path.basename(settings.IN_FILE)))
			xlabel = "Simulated Time Distance"
			title = "Simulation(a=%.3f,e=%.3f,l=%.2f,g=%.2f)"%(settings.INDEL, settings.EPSILON, settings.LAMBDA, settings.GAMMA)
			plotSimulationResult(t, t, oTReal, oIndelibleDif, np.array(eInDel), oInDel, np.array(eMatch), np.array(match), xlabel, fName, title)
	except IOError as e:
                print 'I/0 error{0}:{1}'.format(e.errno, e.strerror)


def findPairBrace(s):
        sLen = len(s)
        right = []
        start = 0
        for i in range(1,sLen):
                while s[i] != '(': break
                if s[i] == '(':
                        right.append(i)
                if len(right) == 1:
                        start = right[0]
                if s[i] == ')':
                        right.pop()
                        if 0 == len(right):
                                return (start,i)
def findTaxaAndLength(s):
        taxa = []
        time = []
        if s[1] == "(":
                index=findPairBrace(s)
                taxa = s[index[0]:index[1]+1]
                time = re.match('^:(\d+\.?\d+),?.+',s[index[1]+1:]).group(1)
        else:
                m = re.match('^\(?(?P<name>[^:]+):(?P<time>\d+\.?\d+).*', s)
                taxa = m.group('name')
                time = m.group('time')
        return (taxa.strip(','), time)
	
def treeSimulation(outFile, seq=None):
	isTaxa = lambda s: (len(re.findall('\(', s)) < 1)
	sTree = settings.TREE
	if None != seq:
		sTree = seq.name
	(left, t1) = findTaxaAndLength(sTree)
	(right, t2) = findTaxaAndLength(sTree[len(left)+len(t1)+2:])
        pair_hmm = pHMM(float(t1), float(t2), left, right)

	rootOrLen = settings.LENGTH
	if None != seq:
		rootOrLen = copy.deepcopy(seq)
		rootOrLen.removeGap()

        sequences = pair_hmm.generateSeq(rootOrLen)
	for i in range(2):
		printFastaFile([sequences[i]], outFile)	
	
	for i in range(2):
		if not isTaxa((left, right)[i]):
			treeSimulation(outFile, sequences[i])	
	return 0

def checkIndelEffect():
	iFile = settings.IN_FILE

	indelArr = [ 0.02 + i * 0.01 for i in range(9)]
	accArr = []
	realTime = settings.TIME
	realIndel = settings.INDEL
	t = [ 0.05 + i*0.01 for i in range(36)]
        for a in indelArr:
                start = time.clock()
                settings.INDEL = a
                settings.IN_FILE = re.sub(r'\.fas$', '_a%.3f.fas'%(settings.INDEL),iFile)
                result = simulationAndAlignmentWithVariantDistance()
		accArr.append(result[:])
	fName = "alignmentOverVariantIndelRateGenerated_t%.2f_a%.2f_e%.2f_l%.2f_g%.2f_len%d_rep%d.png"%(realTime, realIndel, settings.EPSILON, settings.LAMBDA, settings.GAMMA, settings.LENGTH, settings.REPLICATE)
	title = "Alignment Result with variant distance parameters"
	plotAlignmentResult(t, accArr, realTime, fName, title, indelArr)

		
def main(argv):
	#random.seed(1)
	validateArgv(argv)
	#Read command line argument to settings
	readArgv(argv)
	if settings.ALIGNMENT:
		print alignmentWithSpecificDistance(settings.TIME)
	elif settings.PCHECK:
		checkAlignmentParameters()
	elif settings.ACHECK:
		checkIndelEffect()
	elif settings.SIMUL:
		replicateSimulation(settings.IN_FILE)
	elif settings.TCHECK:
		simulationAndAlignmentWithVariantDistance(True)
	elif None == settings.TREE:
		checkAccuracySimulation()
	elif len(re.findall('\(',settings.TREE)) < 2:
		print pairSimulationAndAlignment()
	else:
		try:	
			outFile = open(settings.IN_FILE, 'w')
			treeSimulation(outFile)
		except IOError as e:
			print 'I/O error{0},{1}'.format(e.errno, e.strerror)

if __name__ == "__main__":
	main(sys.argv)
