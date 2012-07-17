#!/usr/bin/python
"""
	The file is used to read the parameters from command line to set the global settings. In the future, we maybe allow using a setting file as an input.
"""
import sys
import settings

def readArgv(argv):

	intArguments = {
		'len': 'LENGTH',
	      'iter' : 'ITERATE',
	      'rep'  : 'REPLICATE'
	}
	floatArguments = {
		'l': 'LAMBDA',
		'e': 'EPSILON',
		'a': 'INDEL',
		'g': 'GAMMA',
		 't' : 'TIME',
	    'step' : 'STEP'
	}
	fileArguments = {
		'in' : 'IN_FILE',
		'out': 'OUT_FILE'
	}
	stringArguments = {
		'tree': 'TREE'
	}
	boolArguments = {
		'acheck': 'ACHECK',
		'pcheck': 'PCHECK' }

	for x in range(1, len(argv)):
		key, value = argv[x].split("=",2)
		rValue = None

		if key in intArguments:
			settings.__dict__[intArguments[key]] = int(value)
		elif key in floatArguments:
			settings.__dict__[floatArguments[key]] = float(value)
		elif key in fileArguments:
			settings.__dict__[fileArguments[key]]  = '%s/%s/%s'%(settings.ROOT, 'data', value)
		elif key in stringArguments:
			settings.__dict__[stringArguments[key]] = value.strip()
		elif key in boolArguments:
			if 0 == int(value):
				settings.__dict__[boolArguments[key]] = False
			else:
				settings.__dict__[boolArguments[key]] = True
				
