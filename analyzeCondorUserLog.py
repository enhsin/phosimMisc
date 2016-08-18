import os, sys
import numpy as np

def str2time(ss):
    s1=ss.split('+')
    dd=float(s1[0])
    s2=s1[1].split(':')
    hh=float(s2[0])
    mm=float(s2[1])
    return dd*24.0+hh+mm/60.0

def getchipID():
    chipID=[]
    fp=open("data/lsst/focalplanelayout.txt").readlines()
    for line in fp:
        lstr=line.split()
        addFlag=0
        if "Group0" in line:
            chipID.append(lstr[0])
    return chipID

def jobCompleted(f):
    if not os.path.exists(f):
        return 0
    b=0
    log=open(f).readlines()
    if len(log)>14:
        if "Normal termination" in log[-10] + log[-14]:
            b=1
    return b

def analyzeLog(fn,objid,visit):
    chipID=getchipID()
    outfile=open(fn+'.log','w')
    i=0
    wtArr=[]
    gtArr=[]
    ctArr=[]
    atArr=[]
    for cid in chipID:
        f='work/logs/log_'+objid+'_'+cid+'_E'+visit+'.log'
        if not jobCompleted(f):
            continue
        #print f
        log=os.popen('condor_userlog -total '+f).readlines()
        log1=log[-1].split()
        wallTime=str2time(log1[1])
        goodTime=str2time(log1[2])
        cpuTime=str2time(log1[3])
        allocTime=str2time(log1[4])
        outfile.write('%8s %8.2f %8.2f %8.2f %8.2f\n' % (cid, wallTime, goodTime, cpuTime, allocTime))
        wtArr.append(wallTime)
        gtArr.append(goodTime)
        ctArr.append(cpuTime)
        atArr.append(allocTime)
        i+=1
    outfile.close()
    wtArr=np.array(wtArr)
    gtArr=np.array(gtArr)
    ctArr=np.array(ctArr)
    atArr=np.array(atArr)
    print 'Summery of '+ str(i) + ' files\n'
    print 'Average:'
    print '  Wall Time  Good Time  CPU Usage  Avg Alloc'
    print '%10.1f %10.1f %10.1f %10.1f' % (wtArr.mean(), gtArr.mean(), ctArr.mean(), atArr.mean())
    print
    print 'Min:'
    print '  Wall Time  Good Time  CPU Usage  Avg Alloc'
    print '%10.1f %10.1f %10.1f %10.1f' % (wtArr.min(), gtArr.min(), ctArr.min(), atArr.min())
    print
    print 'Max:'
    print '  Wall Time  Good Time  CPU Usage  Avg Alloc'
    print '%10.1f %10.1f %10.1f %10.1f' % (wtArr.max(), gtArr.max(), ctArr.max(), atArr.max())


objid=sys.argv[1]
fn=sys.argv[2]
visit=sys.argv[3]

analyzeLog(fn,objid,visit)
