import numpy as np
from matplotlib import pyplot

def plotDistance(t, eTReal, oTReal, oIndelibleDif=None, oPhylosimDif=None):
	#colorCycle = ['black', 'red', 'blue']
	#pyplot.gca().set_color_cycle(colorCycle)
	if isinstance(eTReal, np.ndarray):
		pyplot.plot(t, eTReal)
	else:
		pyplot.axhline(y=eTReal, color='r')
	pyplot.plot(t, oTReal)
	pyplot.ylabel("Observed Time Distance")
	legendArr = ['Expected', 'Simulated', 'Indelible']
	if None == oIndelibleDif:
		legendArr = legendArr[0:2]
	else:
		pyplot.plot(t, oIndelibleDif)
	pyplot.legend(legendArr, loc='upper left')

def plotIndelEvents(t, eInDel, oInDel):
	colorCycle = ['black', 'red']
	if isinstance(eInDel, np.ndarray):
		pyplot.plot(t, eInDel)
	else:
		pyplot.axhline(y=eInDel, color='r')
	pyplot.plot(t, oInDel)
	pyplot.ylabel("The number of Indel Events")
	pyplot.legend(['Expected', 'Simulated'], loc="upper left")

def plotMatchEvents(t, eMatch, oMatch):
	if isinstance(eMatch, np.ndarray):
		pyplot.plot(t, eMatch)
	else:
		pyplot.axhline(y=eMatch, color="r")
	pyplot.plot(t, oMatch)
	pyplot.ylabel("The number of Match Events")
        pyplot.legend(['Expected', 'Simulated'], loc="upper right")

def plotSimulationResult(t, eTReal, oTReal, oIndelibleDif, eInDel, oInDel, eMatch, oMatch, xlabel, fileName, title):
	pyplot.figure(1,(10,10))
	pyplot.title(title)
        pyplot.subplot(311)
        plotDistance(t, eTReal, oTReal, oIndelibleDif)
        pyplot.subplot(312)
        plotIndelEvents(t, eInDel, oInDel)
        pyplot.subplot(313)
        plotMatchEvents(t, eMatch, oMatch)
	pyplot.xlabel(xlabel)
        pyplot.savefig('./image/'+fileName, dpi=300)
	pyplot.close()
        #pyplot.show()

def plotAlignmentResult(t, accuracy,realTime, fileName, title, indel=None):
	pyplot.figure(2,(14,10))
	pyplot.gca().set_color_cycle(['black','red','blue','yellow','green','purple', 'k', 'c', 'm'])
	for a in accuracy:
        	pyplot.plot(np.array(t), np.array(a))
		pyplot.plot(realTime, a[t.index(realTime)],'ro')
	if None != indel:
		pyplot.legend(indel, loc='center left', bbox_to_anchor=(1, 0.5))
        #pyplot.plot(np.array(t), np.array(accuracy))
	pyplot.title(title)
        pyplot.xlabel("distance")
        pyplot.ylabel("accuracy")
        pyplot.xlim(t[0], t[len(t)-1])
        pyplot.ylim(np.array(accuracy).min()-0.01, 1)
	'''
        pyplot.ylim(np.array(accuracy).min()-0.1, 1)
        maxAcc = np.array(accuracy).max()
        pyplot.annotate('local max(%.2f, %.4f)'%(t[accuracy.index(maxAcc)], maxAcc), xy=(t[accuracy.index(maxAcc)], maxAcc), xycoords='data', xytext=(0.6, 0.9), 
textcoords='axes fraction', arrowprops=dict(facecolor="red", shrink=0.02))
        pyplot.annotate('real(%.2f,%.4f)'%(realTime, accuracy[t.index(realTime)]), xy=(realTime, accuracy[t.index(realTime)]), xycoords='data', xytext=(0.5, 0.5)
, textcoords='axes fraction', arrowprops=dict(facecolor="black", shrink=0.02))
	'''
        pyplot.savefig('./image/'+fileName, dpi=200)
	pyplot.close()
        #pyplot.show()
 
