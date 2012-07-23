import os, sys, random, commands, Queue, time, optparse
from pair_hmm import *
from threading import ThreadError
from threading import Thread
from get_free_nodes import get_free_nodes

# parse command lien arguments
'''
parser = optparse.OptionParser(version = "%prog 1.1 on 23.07.2012 by Chunxiang Li", description = "python script for multiple thread pair-wise alignment")
parser.add_option("-i", "--infname", dest="IN_FILE", default="example.fas", help="result file name for simulation, and input file name for alignment at the same time")
parser.add.option("-a","--indelRate", dest="infile")
'''
#global parameters
workingDir = os.getcwd()
job_queue = Queue.Queue()
class Worker(Thread):
	def __init__(self, job_queue, node):
		Thread.__init__(self)
		self.job_queue = job_queue
		self,node = node
	def run(self):
		all_done = 0
		while not all_done:
			try:
				job = self.job_queue.get(0)
				time.sleep(random.randint(5000,6000)/1000.0))
				single_thread(self.node, job)
			except Queue.Empty:
				all_done = 1
	pass

def single_thread(node, job):
	(job_id, t) = job
	print("%s->%s"%(node, job_id))
	try:
		accuracy = commands.getoutput("ssh %s.hpc.cs.helsinki.fi './%s/pair_hmm.py len=%s t=%.2f a=%s e=%s g=%s l=%s rep=%s tree=%s in=%s'"%(node, workingDir, options.length, options.indel, t, options.extension, options.gamma, options.Lambda, options.rep, options.tree, options.ifname))
		os.system("echo '-->realignment job %d end on node %s >>%s'"%(job_id, node, os.path.splitext(options.ifname)[0], logFile))
		os.system("echo '%.4f' >> result_for_%s_with_t%.2f"%(accuracy, os.path.splitext(options.ifname)[0], t))
	except:
		job_queue.put((job))
	time.sleep(random.randint(1000, 5000)/1000.00))
	pass

def simulationAndRealignment():
	#simulation
	os.system("./%s/pair_hmm.py len=%s t=%s a=%s e=%s g=%s l=%s rep=%s tree=%s in=%s >> %s"%(workdingDir, options.length, options.time, options.indel, options.extension, options.gamma, options.Lambda, options.rep, options.tree, options.ifname, logFile))
	
	#realignment
	#get job queue
	for i in range(36):
		job_queue.put(i, 0.05+0.01*i)
	os.system("echo 'processing"%(logFile))
	while not job_queue.empty():
		os.system("echo '---->%d%s >> %s"%(i,cluster[i%len(cluster)], logFile)) 
		t = Worker(job_queue, cluster[i%len(cluster)])
		time.sleep(1)
		try:
			t.start()
			threads.append(t)
		except ThreadError:
			os.system("echo '\t\tError: thread error caught!'>>%s"%(logFile))
	for t in threads:
		t.join()
	#combine results
	os.system("echo 'combining results>>%s"%(logFile))
	for i in range(36):
		jobFName = "result_for_%s_with_t%.2f"%(os.path.splitext(options.ifname)[0], 0.05+0.01*i)
		accuracy.append(real(commands.getoutput("less data/result")))

if __name__ == "__main__":
	global logFile
	validateArgv(sys.argv)
	readArgv(sys.rgv)
	#print get_free_nodes()


