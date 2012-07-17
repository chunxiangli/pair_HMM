#!/usr/bin/python
#from model import *
#from pair_HMM import *
#children = [seqDNA('x'),seqDNA('y')]
#newChildren = [seqDNA('x\''),seqDNA('y\'')]

#children.insertSeq("T-GG")
#children.insertSeq("-C--")

#newChildren.insertSeq("TGG",[1,3,4])
#newChildren.insertSeq("C--",[2])
'''
import re
sTree = "((t1:0.1,(t6:0.1,t2:02):0.2):0.1,(t3:0.2,t4:0.1):0.2)"
sTree = "(t1:0.1,(t6:0.1,t2:02):0.2)"
sTree = "(t1:0.2,t3:0.2)"
sTree = "(((t1:0.1,t7:0.1):0.2,(t6:0.1,t2:02):0.2):0.1,(t3:0.2,t4:0.1):0.2)"
sTree = "(t1:0.2,t3:0.2)"

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
			
		

if __name__ == "__main__":
	#taxa = findPairBrace(sTree)
	#print sTree[taxa[0]:taxa[1]+1]
	#print sTree[taxa[1]+1:]
	#print re.match('^:(\d+\.?\d+),?.+',sTree[taxa[1]+1:]).group(1)
	(taxa, time) = findTaxaAndLength(sTree)
	print (taxa, time)
	print sTree[len(taxa)+len(time)+2:]
	print findTaxaAndLength(sTree[len(taxa)+len(time)+2:])
	print findTaxaAndLength(taxa)
'''
import os
class basic:
	a = 1
	def __init__(self,t):
		a = t
	@staticmethod
	def bPrint(t):
		print t
	def testPrint(self,t):
		self.bPrint(t)
	def getStateDistribution(self):
		transMatrix=[[0.1,0.2,0.7],[0.3,0.4,0.3],[0.6,0.4,0]]
                f = open('input', 'w')
		f.writelines('\n'.join( ' '.join(str(j) for j in i) for i in transMatrix))
                #temp = ' '.join((str(j) for j in i) + '\n' for i in self.transMatrix)
                #f.write(' '.join(str(j) for j in i) + '\n' for i in self.transMatrix)
                f.close()
                os.system("./stateDistribution.o >result")
               	f = open('result', 'r')
                s = f.read()
                f.close()
                return [float(i) for i in s[0].strip('\n').split()]


if __name__ == "__main__":
	for j in range(5):
		t  = basic(j)
		for i in range(10):
			t.getStateDistribution()
		
