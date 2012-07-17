#!/usr/bin/python
import sys
import random

MINFLOAT = 0.000001

def randomIndex(freqArr):
	f = random.random()
	s = 0.0
	iFLen = len(freqArr)
	for i in range(iFLen):
		s += freqArr[i]
		if (s - f) > MINFLOAT:
			return i
	return random.randint(0,iFLen-1)
