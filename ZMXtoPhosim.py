import sys, os, commands
import numpy as np
from decimal import *

class Surface(object):
    def __init__(self,curv,disz,z,rout,rin,coni,an):
        self.curv=curv
        self.disz=disz
        self.z=z
        self.rout=rout
        self.rin=rin
        self.coni=coni
        self.an=an

def findSurface(zmx,surf,surf0,flt,zprev=0.0):
    getcontext().prec = 20
    an=np.zeros(16)
    coni=0.0
    curv=0.0
    disz=[Decimal("0") for i in range(100)]
    pzup=[[Decimal("-1"),Decimal("-1"),Decimal("-1")] for i in range(100)]
    rout=0.0
    rin=0.0
    readCurv=False
    #read DISZ, update THIC, CRVT
    for line in open(zmx).readlines():
        if 'SURF ' in line:
            s=int(float(line.split()[1]))
        elif 'DISZ' in line:
            if line.split()[1]!='INFINITY':
                #disz[s]=float(line.split()[1])
                disz[s]=Decimal(line.split()[1])
        elif 'THIC' in line:
            l=line.split()
            if float(l[2])-1==flt:
                #disz[int(l[1])]=float(l[3])
                disz[int(l[1])]=Decimal(l[3])
        elif 'CRVT' in line:
            l=line.split()
            if float(l[2])-1==flt and float(l[1])==surf:
                curv=-1/float(l[3])
                readCurv=True
        elif 'PZUP' in line:
            l=line.split()
            #pzup[s,0]=float(l[1])
            #pzup[s,1]=float(l[2])
            #pzup[s,2]=float(l[3])
            pzup[s][0]=Decimal(l[1])
            pzup[s][1]=Decimal(l[2])
            pzup[s][2]=Decimal(l[3])

    #PZUP
    for i in range(surf):
        #if pzup[i,0]>0:
            #disz[i]=disz[int(pzup[i,0])]*pzup[i,1]+pzup[i,2]
        if pzup[i][0]>0:
            disz[i]=disz[int(pzup[i][0])]*pzup[i][1]+pzup[i][2]

    z=Decimal("0.0")
    if surf>surf0:
        for i in range(surf):
            if i==surf0:
                z0=z
            z-=disz[i]
        z-=z0

    #read asphere
    found=False
    for line in open(zmx).readlines():
        if 'SURF '+str(surf) in line:
            found=True
            continue
        if found:
            if 'CURV' in line and readCurv==False:
                if float(line.split()[1])!=0:
                    curv=-1/float(line.split()[1])
            elif 'CLAP' in line or 'FLAP' in line:
                rout=float(line.split()[2])
                rin=float(line.split()[1])
            elif 'DIAM' in line:
                rout=float(line.split()[1])
                #print surf, rout
            elif 'CONI' in line:
                coni=float(line.split()[1])
            elif 'PARM' in line:  #parm 1: a2, parm 2: a4
                lstr=line.split()
                an[int(float(lstr[1]))*2-1]=float(lstr[2])/1e3
            elif 'SURF' in line or 'CONF' in line:
                    break

    surface=Surface(curv,z-zprev,z,rout,rin,coni,an)
    return surface


def printOptics(output,surface,name,typ,coating,medium,flt):
    out=open(output,'a')
    if typ=='det' and surface.rout==0:
        surface.rout=400.0
    print '%d %10s %30.20f %30.20f' % (flt,name,surface.disz, surface.z)
    out.write('%-7s %7s %12.5e %20.13e %12.5e %12.4e %12.4e '  % (name,typ,surface.curv,surface.disz,surface.rout,surface.rin,surface.coni))
    for i in range(8):
        out.write('%12.4e ' % (surface.an[i+2]))
    if coating!='none' and typ=='filter':
        coating='filter_'+str(flt)+'.txt'
    out.write('%s %s\n' % (coating,medium))
    out.close()


def checkOptics():
    lsstDir='/home/abby/PROGRAM/phosim/data/lsst/'
    colName=['Name','Type','Curvature','dz','Rout','Rin','Conic','a3','a4','a5','a6','a7','a8','a9','a10','Coating','Medium']
    for i in range(6):
        print 'Differences in optics_'+str(i)+'.txt'
        data0=[]
        for line in open(lsstDir+'optics_'+str(i)+'.txt').readlines():
            if '#' in line:
                continue
            data0.append(line)

        k=0
        for line in open('optics_'+str(i)+'.txt').readlines():
            l=line.split()
            l0=data0[k].split()
            for j in range(len(l)):
                diff=False
                if j in [0,1,15,16]:
                    if l[j]!=l0[j]: diff=True
                else:
                    if float(l0[j])==0:
                        if np.abs(float(l[j])-float(l0[j]))>1e-3: diff=True
                    else:
                        if np.abs((float(l[j])-float(l0[j]))/float(l0[j]))>1e-3: diff=True
                if diff:
                    print '%10s %10s %10s %10s' % (l0[0], colName[j], l[j], l0[j])
            k+=1
        print ''



zmxFile=sys.argv[1]
inputList=sys.argv[2]
fltNum=int(float(commands.getoutput('grep MNUM '+zmxFile+' | awk \'{print $2}\'')))

zmxSurface=[]
name=[]
typ=[]
coating=[]
medium=[]
for line in open(inputList).readlines():
    l=line.split()
    zmxSurface.append(int(float(l[0])))
    name.append(l[1])
    typ.append(l[2])
    coating.append(l[3])
    medium.append(l[4])


for flt in range(fltNum):
    output='optics_'+str(flt)+'.txt'
    try:
        os.remove(output)
    except OSError:
        pass

    zprev=Decimal("0.0")
    for i in range(len(name)):
        surface=findSurface(zmxFile,zmxSurface[i],zmxSurface[0],flt,zprev)
        zprev=surface.z
        printOptics(output,surface,name[i],typ[i],coating[i],medium[i],flt)
#checkOptics()
