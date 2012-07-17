#!/usr/bin/python
import os
import sys
import settings
import re
import random
import copy
import math
import time
from input import *
from model import *
from utils.myprint import *
from utils.tools import *
from utils.myplot import *
import numpy as np
from matplotlib import pyplot
#from scipy import stats


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
		site = real[0].site[index]
		if settings.GAP == real[0].seq[index]:
			pIndex = path[1].site.index(site)
			if settings.GAP == path[0].seq[pIndex]:
				rNum += 1
		elif settings.GAP == real[1].seq[index]:
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

	
def pairSimulation(outFile, rep=1, plot=True):
	if None == settings.TREE:
		(left, t1) = ('A', (settings.TIME/2))
		(right, t2) = ('B', (settings.TIME/2))
	else:
		(left, t1) = findTaxaAndLength(settings.TREE)	
		(right, t2) = findTaxaAndLength(settings.TREE[len(left)+len(t1)+2:])
		t1 = float(t1)
		t2 = float(t2)
	pair_hmm = pHMM(t1, t2, left, right)
	
	p = []
	indelArr = []
	matchArr = []
	for i in range(rep):
		children  = pair_hmm.generateSeq(settings.LENGTH)
		printFastaFile(children, outFile)
                if pair_hmm.observedDif >= 0.75:
			print "To large difference"
			p.append(-1)
		else:
			p.append(pair_hmm.observedDif)
                indelArr.append(pair_hmm.indelNum)
                matchArr.append(pair_hmm.matchNum)
	if plot:
		eInDel = pair_hmm.getExpectedInDel()
		eMatch = pair_hmm.getExpectedMatch()
		fName = "estimate_distance.len%d_a%.2f_e%.2f_rep%d.png"%(settings.LENGTH, settings.INDEL, settings.EPSILON, settings.REPLICATE)
		xlabel = "Replicate Sequence"
		title = "Simulation(a=%.3f,e=%.3f,l=%.2f,g=%.2f)"%(settings.INDEL, settings.EPSILON, settings.LAMBDA, settings.GAMMA)
        	plotSimulationResult(np.arange(settings.REPLICATE), t1+t2, d2t(p2d(np.array(p))), None, eInDel, np.array(indelArr), eMatch, np.array(matchArr), xlabel, fName, title)
	
	return (t1 + t2)

def alignment(inFile, outFile, t=0):
	sequences = []
	alignment = []

	s = inFile.readline()
	for i in range(2):
		while '' != s and '>' != s[0]: s = inFile.readline()
		if '' == s:
			print "There isn\'t sequence in the fasta file %s."%(settings.IN_FILE)
			inFile.close()
			sys.exit(0)
		name = s[1:].split('|')[0]
		seq = ''
		s = inFile.readline()
		while '' != s and '>' != s[0]: 
			seq += s
			s = inFile.readline()
		if '' == s and '' == seq:
			print "There is only one sequence in the fasta file %s."%(settings.IN_FILE) 
			inFile.close()
			sys.exit()
		sequences.append(seqDNA(name, re.sub(r'\n', '', seq)))
	inFile.seek(inFile.tell() - len(s))
	if 0 == t:
		t1 = float(re.search(':(\d+\.?\d+)\|?', sequences[0].name).group(1).strip())
		t2 = float(re.search(':(\d+\.?\d+)\|?', sequences[1].name).group(1).strip())
	else:
		t1 = settings.TIME / 2
		t2 = settings.TIME / 2

	pair_hmm = pHMM(t1, t2)
	alignment = copy.deepcopy(sequences)
	#printPairAlign([alignment[0].seq, alignment[1].seq])
	for c in alignment:
		c.removeGap()

	start = time.clock()
	alignment[0],alignment[1] = pair_hmm.alignSeq(alignment)
	end = time.clock()
	#printPairAlign([alignment[0].seq, alignment[1].seq])
	printFastaFile(alignment, outFile)
	accuracy = fidelity(sequences, alignment)
	for i in range(2):
		sequences[i].removeGap()
		alignment[i].removeGap()
		a1 =''.join(sequences[i].seq)
		a = ''.join(alignment[i].seq)
		if a1 != a:
			print '%s\n%s'%(a1,a) 

	observedP = pair_hmm.observedDif
	return (observedP, accuracy)

def pairSimulationAndAlignment(plot=True):
	replicate = settings.REPLICATE 
	P = []
	t = []
	accuracy = []
	realTime = 0
	try:
		realTime = replicateSimulation(settings.IN_FILE, replicate)

		f1 = open(settings.IN_FILE, 'r')
		for j in range(2):
			f2 = open(re.sub(r'\.fas$', 'a%.3f_t%.3f.fas'%(settings.INDEL, settings.TIME), settings.OUT_FILE),'w')
			settings.TIME = 0.05 + j*0.01
			t.append(settings.TIME)
			accArr = []	
			for i in range(replicate):
				(a, b) = alignment(f1, f2, 1)
				#P.append(a)
				accArr.append(b)
			f1.seek(0)
			accuracy.append(np.array(accArr).mean())
			f2.close()
		f1.close()
		if plot:
			fName = "realignment_t0.05:0.4_a%.3f_e%.3f_len%d_rep%d.png"%(settings.INDEL, settings.EPSILON, settings.LENGTH, settings.REPLICATE)
			title = "Alignment(a=%.3f,e=%.3f,l=%.2f,g=%.2f)"%(settings.INDEL, settings.EPSILON, settings.LAMBDA, settings.GAMMA)
			plotAlignmentResult(t, [accuracy], realTime ,fName, title)

	except IOError as e:
        	print 'I/0 error{0}:{1}'.format(e.errno, e.strerror)
	return (t, accuracy, realTime)

def replicateSimulation(outFile, replicate=1, plot=True):
	of = open(outFile, 'w')
	realTime = pairSimulation(of, replicate, plot)
	of.close()	
	return realTime

def checkAlignmentParameters():
	indel = settings.INDEL
	t = settings.TIME
	replicateSimulation(settings.IN_FILE, settings.REPLICATE, False)
	
	indelArr = [ 0.02 + i * 0.01 for i in range(9)]
	tArr = [0.05 + i * 0.01 for i in range(36)]
	accuracy = []

	iF = open(settings.IN_FILE, 'r')
	for i in range(9):
		settings.INDEL = 0.02 + i * 0.01
		accArr = []
		for j in range(36):
			ofName = re.sub(r'\.fas$', 'a%.3f_t%.3f.fas'%(settings.INDEL, settings.TIME), settings.OUT_FILE)
			#of = open(re.sub(r'\.fas$', 'a%.3f_t%.3f.fas'%(settings.INDEL, settings.TIME), settings.OUT_FILE),'w')
			of = open(ofName, 'w')
			settings.TIME = 0.05 + j * 0.01
			acc = []
			for k in range(settings.REPLICATE):
				(a, b) = alignment(iF, of, 1)
				acc.append(b)
			of.close()
			iF.seek(0)
			accArr.append(np.array(acc).mean())
		accuracy.append(accArr)

	fName = "alignmentOverVariantParameters_generate(t%.3f,a%.3f,e%.3f,l%.3f,g%.3f).fas"%(t, indel, settings.EPSILON, settings.LAMBDA, settings.GAMMA)
	title = "Alignment(t=%.3f, a=%.3f,e=%.3f,l=%.2f,g=%.2f)"%(t, indel, settings.EPSILON, settings.LAMBDA, settings.GAMMA)
	plotAlignmentResult(tArr, accuracy, t ,fName, title, indelArr)
		
	iF.close()

def runIndelible(start, step, times, subRate, length):
	observedDif = []
	os.system("./runIndelible %f %f %d %f %d>>log.txt"%(start, step, times, subRate, length))
	for i in range(times):
		iF = "/home/czli/Downloads/INDELibleV1.03/data/output%.3f.fas"%((start + i*step)/2)
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
	localIter = 1 
	try:
		start = settings.TIME
		f1 = open(settings.IN_FILE, 'w')
		for i in range(settings.ITERATE):
			t.append(settings.TIME)
			sumD = 0
			wIter = 0
			indelArr = []
			matchArr = []

			pair_hmm = pHMM(settings.TIME/2, settings.TIME/2, "A", "B")

			for j in range(localIter):
        			children  = pair_hmm.generateSeq(settings.LENGTH)
        			printFastaFile(children, f1)
				p = pair_hmm.observedDif
				if p >= 0.75:
					wIter += 1
				else:
					sumD += p2d(p)	
				indelArr.append(pair_hmm.indelNum)
				matchArr.append(pair_hmm.matchNum)

			eInDel.append(pair_hmm.getExpectedInDel())
			indel.append(np.array(indelArr).mean())
			match.append(np.array(matchArr).mean())
			eMatch.append(pair_hmm.getExpectedMatch())
			d.append(sumD/(localIter - wIter))
			settings.TIME += settings.STEP		
		f1.close()
		indelDif = runIndelible(start, settings.STEP, settings.ITERATE, settings.INDEL, settings.LENGTH)
		#phylosimDif = runPhylosim(start, settings.STEP, settings.ITERATE, settings.INDEL, settings.LENGTH)
      		
		if 1 != settings.ITERATE:
			#print d
			t = np.array(t)
			oDReal = np.array(d)
			oTReal = d2t(oDReal)
			oInDel = np.array(indel)
			oIndelibleDif = d2t(p2d(np.array(indelDif)))
			#oPhylosimDif = d2t(p2d(np.array(phylosimDif)))
			#s = stats.linregress(t, oTReal)
			#lTReal = s[0]*t+s[1]
			fName = "estimate_distance.len%d_a%.2f_e%.2f_iter%d.png"%(settings.LENGTH, settings.INDEL, settings.EPSILON, settings.ITERATE)
			xlabel = "Simulated Time Distance"
			title = "Simulation(a=%.3f,e=%.3f,l=%.2f,g=%.2f)"%(settings.INDEL, settings.EPSILON, settings.LAMBDA, settings.GAMMA)
			plotSimulationResult(t, t, oTReal, oIndelibleDif, np.array(eInDel), oInDel, np.array(eMatch), match, xlabel, fName, title)
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
        oFile = settings.OUT_FILE

	indelArr = []
	accArr = []
	realTime = settings.TIME
	t = []
        for a in range(9):
        	print "a=%d"%(a)
                start = time.clock()
                settings.INDEL = 0.02 + a * 0.01
		indelArr.append(settings.INDEL)
                settings.IN_FILE = re.sub(r'\.fas$', '_a%.3f.fas'%(settings.INDEL),iFile)
                settings.OUT_FILE = re.sub(r'\.fas$', '_a%.3f.fas'%(settings.INDEL),oFile)
                result = pairSimulationAndAlignment(False)
		accArr.append(result[1])
		t = result[0]
	fName = "alignmentOverVariantIndelRate_t%.3f_e%.3f_l%.3f_g%.3f.png"%(realTime, settings.EPSILON, settings.LAMBDA, settings.GAMMA)
	title = "Alignment Result with variant distance parameters"
	plotAlignmentResult(t, accArr, realTime, fName, title, indelArr)

		
def main(argv):
	#random.seed(1)
	validateArgv(argv)
	#Read command line argument to settings
	readArgv(argv)
	if settings.PCHECK:
		checkAlignmentParameters()
	elif False != settings.ACHECK:
		checkIndelEffect()
	elif None == settings.TREE:
		checkAccuracySimulation()
	elif len(re.findall('\(',settings.TREE)) < 2:
		pairSimulationAndAlignment()
	else:
		try:	
			outFile = open(settings.IN_FILE, 'w')
			treeSimulation(outFile)
		except IOError as e:
			print 'I/O error{0},{1}'.format(e.errno, e.strerror)

if __name__ == "__main__":
	main(sys.argv)
