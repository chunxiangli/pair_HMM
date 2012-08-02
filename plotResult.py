#!/usr/bin/python
import optparse
import os
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
(options, arguments) = parser.parse_args()

if not options.dataDir:
	print "Please input the result directory using -d"
	sys.exit(0)

ROOT = os.getcwd()
ifname = "%s/data/%s/finalResult"%(ROOT, options.dataDir)
f = open(ifname, 'r')
s = f.read()
result = [ float(i) for i in s.strip().split()]
f.close()
t = [ 0.05 + 0.01 * i for i in range(56)]
fName = 'alignmentWithVariantDistance_over_%s.png'%(options.dataDir)
title = 'alignment(a=%s, e=%s, l=%s, g=%s, len=%s, rep=%s)'%(options.indel, options.extension, options.gamma, options.Lambda, options.length, options.rep)
plotAlignmentResult(t, [result], float(options.time), fName, title)
#plot probability score
ifname = "%s/data/%s/finalScore"%(ROOT, options.dataDir)
f = open(ifname, 'r')
s = f.read()
result = [float(i) for i in s.strip().split()]
f.close()
fName = 'probability_alignmentWithVariantDistance_over_%s.png'%(options.dataDir)
plotAlignmentResult(t, [result], float(options.time), fName, title, None, 'probability(log)')
