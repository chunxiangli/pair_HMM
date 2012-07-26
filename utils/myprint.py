#!/usr/bin/python
import sys
import re
import settings
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


def printPairAlign(Pair):
	assert 2 == len(Pair),"Not a pair of sequences!"

	for i in Pair:
		assert isinstance(i,str) or isinstance(i,list),"Sequence format error!"
	sLen = len(Pair[0])
	arr = [[]]*3
	arr[2] = ['|']*sLen
	arr[:2] = Pair

	for i in range(len(Pair)):
		if isinstance(Pair[i],list):
			arr[i] = ''.join(Pair[i])

	for i in range(sLen):
		if "-" == Pair[0][i] or "-" == Pair[1][i]:
			arr[2][i] = ' '
	if sLen > 100:
		gen = sLen /100
		for i in range(gen):
			start = i*100
			end = (i+1)*100
			print arr[0][start:end]
			print ''.join(arr[2][start:end])
			print arr[1][start:end]
			print "\n"
		if sLen - gen*100 > 0:
			start = gen*100
			print arr[0][start:]
			print ''.join(arr[2][start:])
			print arr[1][start:]
	else:	
		print arr[0]
		print ''.join(arr[2])
		print arr[1]

'''	
def printFastaFile(Seq, outFile):
	for s in Seq:
		outFile.write('>%s\n'%(s.name))
		sLen = s.length / settings.FASTA_LINE_LENGTH
		start = 0
		if s.length > settings.FASTA_LINE_LENGTH:
			for i in range(sLen):
				end = start +  settings.FASTA_LINE_LENGTH
				outFile.write(''.join(s.seq[start:end])+'\n')
				start = (i + 1) * settings.FASTA_LINE_LENGTH
		outFile.write(''.join(s.seq[start:])+'\n')
'''

def seq2Record(seq):
	return SeqRecord(Seq(''.join(seq.seq)), id=seq.name, name='', description='')

def printFastaFile(seq, outFile):
	if not isinstance(seq[0], list):
		SeqIO.write([ seq2Record(s) for s in seq], outFile, "fasta")
	else:
		for s in seq:
			SeqIO.write([seq2Record(i) for i in s], outFile, 'fasta')

def printPhylipFile(seq, outFile):
	SeqIO.write([ seq2Record(s) for s in seq], outFile, "phylip")

	
if __name__ == "__main__":
	printSeq(['a','b','c'])
