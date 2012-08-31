#!/usr/bin/python
import os
import sys
import random
import commands
import Queue
import time
import optparse
import numpy as np
import multiprocessing
from pair_hmm import *
from threading import ThreadError
from threading import Thread
from get_free_nodes import get_free_nodes 

# parse command lien arguments
parser = optparse.OptionParser(version = "%prog 1.1 on 23.07.2012 by Chunxiang Li", description = "python script for multiple thread pair-wise alignment")
parser.add_option("-i", "--in", dest="ifname", default="example.fas", help="result file name for simulation, and input file name for alignment at the same time")
parser.add_option("-o", "--out", dest="ofname", default="example.fas", help="result directory name")
parser.add_option("-t","--time", dest="time", default="0.4", help="expected number of substitution events per sites")
parser.add_option("-e", "--e", dest="extension", default="0.4", help="gap extension probability")
parser.add_option("-l", "--l", dest="Lambda", default="1", help="substitution rate per unit time")
parser.add_option("-g", "--g", dest="gamma", default="0", help="gap extension probability")
parser.add_option("-s", "--slen", dest="length", default="100", help="result alignment length")
parser.add_option("-r", "--rep", dest="rep", default="10", help="replicate num")
parser.add_option("-d", "--dtree", dest="tree", default="(A:0.1,B:0.1)", help="replicate num")
parser.add_option("-m", "--model", dest="model", default="Durbin", help="alignment model name")
parser.add_option("-f", "--f", dest="format", default="fasta", help="input file format")
(options, arguments) = parser.parse_args()
max_thread = 5#int(options.max_thread)
cpu_num = multiprocessing.cpu_count()
seq_num = int(options.rep)
max_node = 130 

block = seq_num/(max_thread*cpu_num)
if block * max_thread * cpu_num < seq_num:
	block += 1

#global parameters
workingDir = os.getcwd()
job_queue = Queue.Queue()
class Worker(Thread):
	def __init__(self, job_queue, node):
		Thread.__init__(self)
		self.job_queue = job_queue
		self.node = node
	def run(self):
		all_done = 0
		while not all_done:
			try:
				job = self.job_queue.get(0)
				time.sleep(random.randint(5000,6000)/1000.0)
				single_thread(self.node, job)
			except Queue.Empty:
				all_done = 1
	pass

def single_thread(node, job):
	(job_id, a, block_id) = job
	print("%s->%s %.2f %d"%(node, job_id, a, block_id))
	try:
		accuracy = commands.getoutput('''ssh -o StrictHostKeyChecking=no %s 'cd Documents/IB/pair_HMM/;./pair_hmm.py len=%s t=%s a=%.2f e=%s g=%s l=%s rep=%s id=%d in=%s out=%s m=%s itype=%s alignment=1 acheck=1' '''%(node,  options.length, options.time, a, options.extension, options.gamma, options.Lambda, options.rep, block_id, options.ifname, options.ofname, options.model, options.format))
		os.system("echo '-->realignment job %d_%d end on node %s' >>%s"%(job_id, block_id,  node, logFile))
		os.system("echo '%s' > %s/result_for_%s_%d_with_a%.2f"%(accuracy, dataDir, os.path.splitext(os.path.basename(options.ifname))[0], block_id, a))
	except:
		job_queue.put((job))
	time.sleep(random.randint(1000, 5000)/1000.00)
	pass

def simulationAndRealignment():
	'''
	#simulation
	print "start simulation"
	os.system("echo 'start simulation' > %s"%(logFile))
	block_size = cpu_num*max_thread
	for block_id in range(block):
		if (block_id + 1) * block_size > seq_num:
			block_size = seq_num - block_id * block_size 
		os.system("./pair_hmm.py len=%s a=%s e=%s g=%s l=%s rep=%d tree='%s' id=%d in=%s simul=1>> %s"%(options.length, options.indel, options.extension, options.gamma, options.Lambda, block_size, options.tree, block_id, options.ifname, logFile))
	'''
	#realignment
	#get job queue
	a_range = 12 
	print "get the job queue"
	job_id = 0
	for i in range(a_range):
		for j in range(block):
			job_queue.put((job_id, 0.05 + 0.05 * i, j))
			job_id += 1
	#processing
	job_size = job_queue.qsize()
        print "processing %d jobs" % (job_size)
	start = time.time()
	os.system("echo 'processing %d jobs' >> %s"%(job_size, logFile))
	threads = []
	cluster = get_free_nodes()[0]
	#cluster=['ukko003.hpc','ukko004.hpc','ukko005.hpc','ukko006.hpc','ukko007.hpc','ukko008.hpc']
	load = 1
	for i in range(load):
		if job_queue.empty():
			break
		nodes_num = 0
		for j in range(max_node):
			if job_queue.empty():
				break
			os.system("echo '---->%d  %s'>> %s"%(j, cluster[j%len(cluster)], logFile)) 
			t = Worker(job_queue, cluster[j%len(cluster)])
			nodes_num += 1
			time.sleep(1)
			try:
				t.start()
				threads.append(t)
			except ThreadError:
				os.system("echo '\t\tError: thread error caught!'>>%s"%(logFile))
		time.sleep(random.randint(5000,6000)/1000.0)
		
		os.system("echo 'using %d nodes in cluster' >>%s"%(nodes_num,  logFile))
	for t in threads:
		t.join()
	#combine results
	os.system("echo 'combining results' >>%s"%(logFile))
	accuracy = []
	score = []
	for i in range(a_range):
		acc = []
		lScore = []
		for j in range(block):
			jobFName = "%s/result_for_%s_%d_with_a%.2f"%(dataDir, os.path.splitext(os.path.basename(options.ifname))[0], j, 0.05+0.05*i)
			s = commands.getoutput("tail -1 %s |tr -d '(|)|'|tr ',' ' '"%(jobFName))
			s = s.split()
			acc.append(float(s[0]))
			lScore.append(float(s[1]))
			#os.system("rm "+jobFName)
		accuracy.append(np.array(acc).mean())
		score.append(np.array(lScore).mean())
	os.system("echo %s > %s/finalResult"%(' '.join([str(i) for i in accuracy]), dataDir))
	os.system("echo %s > %s/finalScore"%(' '.join([str(i) for i in score]), dataDir))
	os.system("echo '%.4gs' >> %s"%(time.time()-start, logFile))

if __name__ == "__main__":
	global dataDir 
	global logFile 
	dataDir = 'data/%s'%(options.ofname)
	if not os.path.exists(dataDir):
		os.mkdir(dataDir)
	logFile = dataDir+'/log.txt'
	simulationAndRealignment()
	#print get_free_nodes()
