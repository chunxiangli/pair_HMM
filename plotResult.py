#!/usr/bin/python
import optparse
import os
import re
import commands
import numpy as np
from utils.myplot import *
parser = optparse.OptionParser(version = "%prog 0.1 on 25.07.2012 by Chunxiang Li", description="python script for plotting simulation and alignment result")
parser.add_option("-d", "--data", dest="dataDir",help="the result directory in data/")
parser.add_option("-a", "--a", dest="indel", default="0.1", help="insertion deletion rate")
parser.add_option("-e", "--e", dest="extension", default="0.4", help="gap extension probability")
parser.add_option("-l", "--lambda", dest="Lambda", default="1", help="substitution rate per unit time")
parser.add_option("-g", "--gamma", dest="gamma", default="0", help="match extension probability")
parser.add_option("-s", "--slen", dest="length", default="1000", help="alignment length")
parser.add_option("-r", "--rep", dest="rep", default="1000", help="replicates number")
parser.add_option("-t", "--time", dest="time", default="0.1", help="simulated time distance")
parser.add_option("-y", "--type", dest="ctype", default="t", help="check time effect or indel effect")
(options, arguments) = parser.parse_args()
def getStd():
	filelist = commands.getoutput('ls ./data/%s|grep result_for'%options.dataDir).split()
	prefix = re.sub('_\d+_','',filelist[0]).split('with_')[0]
	y = [0.05 + 0.05*i for i in range(12)]
	block = int(options.rep)/80
	if block*80 < int(options.rep):
		block += 1
	sStd = []
	aStd = []
	for i in y:
		acc = []
		score = []
		for j in range(block):
			f = open("./data/%s/%s_%d_with_%s%.2f"%(options.dataDir, prefix, j, options.ctype, i), 'r')
			s = re.sub(',', '', f.readlines()[0]).strip('()\n').split()
			acc.append(float(s[0]))
			score.append(float(s[1]))
			f.close()
		aStd.append(np.array(acc).std())	
		sStd.append(np.array(score).std())
	return (aStd, sStd)

if not options.dataDir:
	print "Please input the result directory using -d"
	sys.exit(0)

ROOT = os.getcwd()
ifname = "%s/data/%s/finalResult"%(ROOT, options.dataDir)
f = open(ifname, 'r')
s = f.read()
result = [ float(i) for i in s.strip().split()]
f.close()
t = [ 0.05 + 0.05 * i for i in range(12)]
realT = float(options.time)
if 'a' == options.ctype:
	realT = float(options.indel)
(aStd, sStd) = getStd()
fName = 'alignmentWithVariantDistance_over_%s.png'%(options.dataDir)
title = 'alignment(a=%s, t=%s, e=%s, l=%s, g=%s, len=%s, rep=%s)'%(options.indel, options.time, options.extension, options.Lambda, options.gamma, options.length, options.rep)
plotAlignmentResult(t, [result], realT, fName, title, None, 'accuracy', aStd, options.ctype)
#plot probability score
ifname = "%s/data/%s/finalScore"%(ROOT, options.dataDir)
f = open(ifname, 'r')
s = f.read()
result = [float(i) for i in s.strip().split()]
f.close()
fName = 'probability_alignmentWithVariantDistance_over_%s.png'%(options.dataDir)
plotAlignmentResult(t, [result], realT, fName, title, None, 'probability(log)', sStd, options.ctype)

