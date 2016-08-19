import os, sys
from collections import defaultdict
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

def getTime(log):
    log1=log[-1].split()
    wallTime=str2time(log1[1])
    goodTime=str2time(log1[2])
    cpuTime=str2time(log1[3])
    allocTime=str2time(log1[4])
    return [wallTime, goodTime, cpuTime, allocTime]

def addHost(hostDict,log):
    for l in log[3:]:
        if l=='\n':
            break
        lstr=l.split()
        h=lstr[0]
        t=str2time(lstr[3])
        hostDict[h]+=t

def getHostName(h):
    log=os.popen('nslookup '+h).readlines()
    for l in log:
        if 'name =' in l:
            break
    lstr=l.split()[3].split('.')
    if lstr[2]=='purdue':
        if lstr[1]=='rcac':
            hostName=lstr[0].split('-')[0]
        elif lstr[1]=='bio':
            hostName=lstr[1]
        else:
            hostName=lstr[0]+'.'+lstr[1]
    elif lstr[2]=='wisc' and lstr[1]=='chtc':
        hostName=lstr[1]+'.'+lstr[2]
    else:
        hostName=l
    return hostName



def analyzeLog(fn,objid,visit):
    chipID=getchipID()
    out1=open(fn+'_time.log','w')
    out2=open(fn+'_host.log','w')
    out3=open(fn+'_sum.log','w')
    tArr=[]
    hostDict=defaultdict(float)
    i=0
    for cid in chipID:
        f='work/logs/log_'+objid+'_'+cid+'_E'+visit+'.log'
        if not jobCompleted(f):
            continue
        #print f
        log=os.popen('condor_userlog -total '+f).readlines()
        t1=getTime(log)
        tArr.append(t1)
        out1.write('%8s %8.2f %8.2f %8.2f %8.2f\n' % (cid, t1[0], t1[1], t1[2], t1[3]))
        addHost(hostDict,log)
        i+=1
    out1.close()
    tArr=np.array(tArr)
    out3.write('Summery of '+ str(i) + ' files\n\n')
    out3.write('Average:\n')
    out3.write('  Wall Time  Good Time  CPU Usage  Avg Alloc\n')
    out3.write('%10.1f %10.1f %10.1f %10.1f\n\n' % (tArr[:,0].mean(), tArr[:,1].mean(), tArr[:,2].mean(), tArr[:,3].mean()))
    out3.write('Min:\n')
    out3.write('  Wall Time  Good Time  CPU Usage  Avg Alloc\n')
    out3.write('%10.1f %10.1f %10.1f %10.1f\n\n' % (tArr[:,0].min(), tArr[:,1].min(), tArr[:,2].min(), tArr[:,3].min()))
    out3.write('Max:\n')
    out3.write('  Wall Time  Good Time  CPU Usage  Avg Alloc\n')
    out3.write('%10.1f %10.1f %10.1f %10.1f\n\n' % (tArr[:,0].max(), tArr[:,1].max(), tArr[:,2].max(), tArr[:,3].max()))

    hostSumDict=defaultdict(float)
    hostCountDict=defaultdict(int)
    for h in sorted(hostDict, key=hostDict.get, reverse=True):
        out2.write("%16s %10.1f\n" % (h,hostDict[h]))
        if hostDict[h]>0:
            hostName=getHostName(h)
            hostSumDict[hostName]+=hostDict[h]
            hostCountDict[hostName]+=1
    out2.close()

    out3.write('\n')
    for h in sorted(hostSumDict, key=hostSumDict.get, reverse=True):
        out3.write('%20s %10.1f %3d\n' % (h,hostSumDict[h],hostCountDict[h]))
    out3.close()
    os.system('cat '+fn+'_sum.log')



objid=sys.argv[1]
fn=sys.argv[2]
visit=sys.argv[3]

analyzeLog(fn,objid,visit)
