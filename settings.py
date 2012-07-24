'''
	Config file
'''
import os
ROOT = os.getcwd()
######### Algorithm setting part ####################
# base frequency, which can be modified when run
BASE_FREQUENCY = [1.0/4] * 4

# base pool
BASE = ['T','C','A','G']

# sequence length
LENGTH = 20

#default model JC69
LAMBDA = 0.1

MODEL = "JC69"

'''
	parameters of transition matrix
'''
#indel rate al*t < 0.405, otherwise, the probability of match even smaller than opening a gap
INDEL = 0.8

#match extension probability
GAMMA = 0.5

#gap extension probability
EPSILON = 0.4


#gap wildcard
GAP = "-"

# tree architecture description
#TREE = '(t1:0.1, t2:0.1)'
TREE = None

##########

#evolution distance between pair sequences
TIME = 0

STEP = 0

ITERATE = 1

REPLICATE = 1

ACHECK = False

PCHECK = False

TCHECK = False

ALIGNMENT = False
# output file name
OUT_FILE = ROOT + "/data/realign.fas"

# input file name
IN_FILE = ROOT + "/data/example.fas"

IndelibleDataDir = "/home/czli/Downloads/INDELibleV1.03/data"

# the line length in fasta format
FASTA_LINE_LENGTH = 80
