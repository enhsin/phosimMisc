import sys,os,commands
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from mpl_toolkits.mplot3d import Axes3D

from matplotlib.backends.backend_pdf import PdfPages


class photon(object):
    def  __init__(self):
        self.wavelength = 0
        self.numLayer = 0 
        self.xdir = 0
        self.ydir = 0
        self.time = 0
        self.listX = []
        self.listY = []
        self.listZ = []
        self.listLayer = []
        self.opticX = []
        self.opticY = []
        self.opticZ = []
        self.opticn = []
        self.opticL = []
                
    def printPhoton(self):
        #print "****************************PHOTON*********************************"
        #print "wavelength " , self.wavelength 
        #print "X-direction ", self.xdir
        #print "Y-direction ", self.ydir
        #print "time ", self.time
        #print "number of Layers ", self.numLayer
        #print "*******"  
        #print "Layer\tX\tY\tZ"
        #for i in range(0, self.numLayer):
        #    print self.listLayer[i], self.listX[i], self.listY[i], self.listZ[i]
        #print "******************************************************************"       
        print self.listLayer[8], self.listX[8], self.listY[8], self.listZ[8], self.listZ[7], self.listY[19]

def readEvents(eventFits,verbose=1,addIntercept=0,addOP=0):
    plotEvery=False
    f= fits.open(eventFits)
    event=f[1].data
    numLine = len(event.field(3))
    print numLine
    numLine20=int(numLine/500)
    numLine20=1
    photonList = []
    pre_p = photon()
    landed=0
    j=1
    for i in range(0, numLine):
        #if verbose>0:
        #    if i%numLine20==0:
        #        print j,'/20'
        #        j+=1
        if event.field(3)[i]==0:
            if i!=0:
                pre_p.numLayer=numLayer
                if plotEvery:
                    landed=1
                #pre_p.resize_layer()
                if landed==1:
                    photonList.append(pre_p)
                    if verbose>0:
                        if i%numLine20==0:
                            pre_p.printPhoton()
                    landed=0
            numLayer = 0
            p = photon()
            p.wavelength = event.field(2)[i]
            p.xdir = event.field(0)[i]
            p.ydir = event.field(1)[i]
        elif event.field(3)[i]==1:
            p.time = event.field(0)[i]
        else:
            if np.isnan(event.field(1)[i]):
                p = photon()
                continue
            else:
                if event.field(3)[i]!=303:
                    numLayer += 1
                else:
                    landed=1
                if event.field(3)[i]==114 or event.field(3)[i]==10 or event.field(3)[i]==11:
                    numLayer-=1
                    continue
                if event.field(3)[i]>=400:
                    numLayer-=1
                    if addOP:
                        p.opticZ.append(event.field(0)[i])
                        p.opticn.append(event.field(1)[i])
                        p.opticL.append(event.field(2)[i])
                        if addIntercept:
                            if event.field(3)[i]==411:
                                p.opticZ.append(np.sqrt((p.listZ[-1]-p.listZ[-2])**2.0+(p.listY[-1]-p.listY[-2])**2.0))
                                p.opticn.append(1.0)
                                p.opticL.append(0.0)
                    continue
                p.listX.append(event.field(0)[i])
                p.listY.append(event.field(1)[i])
                p.listZ.append(event.field(2)[i])
                p.listLayer.append(event.field(3)[i])
                if addOP and (event.field(3)[i]<200 or event.field(3)[i]==303):
                    p.opticZ.append(0)
                    p.opticn.append(0)
                    p.opticL.append(0)
                if addIntercept:
                    if event.field(3)[i]==211:  #detector
                        dy=p.listY[-1]-p.listY[-2]
                        dz=p.listZ[-1]-p.listZ[-2]
                        dx=p.listX[-1]-p.listX[-2]
                        v=-p.listY[-1]/dy
                        p.listX.append(p.listX[-1]+v*dx)
                        p.listY.append(0.0)
                        p.listZ.append(p.listZ[-1]+v*dz)
                        p.listLayer.append(213)

                pre_p = p
                if i==numLine-1:
                    pre_p.numLayer=numLayer
                    if plotEvery:
                        landed=1
                    if landed==1:
                        photonList.append(pre_p)
                        if verbose>0:
                            if i%numLine20==0:
                                pre_p.printPhoton()
                        landed=0
    print len(photonList)
    return photonList

def findHeight(opticsFile,device):
    height=0.0
    Z=0.0
    for line in open(opticsFile).readlines():
        if "#" in line:
            continue

        lstr=line.split()
        a=np.zeros(13)
        for i in range(13):
            a[i]=float(lstr[i+2])
        height+=a[1]

        if device in line:
            radiusofcurv=-a[0]
            outerRadius=a[2]
            innerRadius=a[3]
            conic=a[4]
            R=outerRadius
            if radiusofcurv==0:
                Z=R-R+height
            else:
                Z=height-R**2/radiusofcurv/(1.0+np.sqrt(1.0-(conic+1.0)*R**2/radiusofcurv**2))
                for i in range(3,11):
                    Z-=a[i+2]*1e3*R**i
            break

    return Z


def plotOptics(opticsFile):
    height=0.0
    for line in open(opticsFile).readlines():
        if "#" in line:
            continue

        lstr=line.split()
        a=np.zeros(13)
        for i in range(13):
            a[i]=float(lstr[i+2])
        height+=a[1]

        surfaceName=lstr[0]
        if surfaceName=='none':
            continue
        radiusofcurv=-a[0]
        outerRadius=a[2]
        innerRadius=a[3]
        innerRadius=0.0
        conic=a[4]
        #u = np.linspace(0, 2 * np.pi, 30)
        u = np.array([np.pi/2, np.pi*3/2.])
        v = np.linspace(innerRadius, outerRadius, 20)
        X = np.outer(v,np.cos(u))
        Y = np.outer(v,np.sin(u))
        R=np.sqrt(X**2 + Y**2)
        if radiusofcurv==0:
            Z=R-R+height
        else:
            Z=height-R**2/radiusofcurv/(1.0+np.sqrt(1.0-(conic+1.0)*R**2/radiusofcurv**2))
            for i in range(3,11):
                Z-=a[i+2]*1e3*R**i
        plt.plot(Z,Y,'-k',linewidth=1.5)
        if surfaceName[-1]=='E':
            plt.text(np.mean(Z),1.15*np.min(Y),surfaceName,fontsize=10)
        else:
            plt.text(np.mean(Z),1.1*np.max(Y),surfaceName,fontsize=10)


    #ep=0.0
    #plt.plot([ep,ep],[-5000,5000],':r')
    #plt.text(ep+40,4000,'entrance pupil',fontsize=10)

    u = np.linspace(-np.pi/4, np.pi/4, 30)
    epR=2738.30878556
    #plt.plot(epR*np.cos(u)+4428.73, epR*np.sin(u),'--k')
    #off-axis
    ag=1.7*np.pi/180.0
    y0=epR*np.sin(u)
    x0=epR*np.cos(u)
    dy=1.7* (180000.0/1000)
    dx=4428.73415441
    xp=np.cos(ag)*x0+np.sin(ag)*y0
    yp=-np.sin(ag)*x0+np.cos(ag)*y0
    plt.plot(xp+dx,yp+dy,'--k')
    plt.plot(x0+dx,y0,'--r')
    plt.text(epR+dx,500,'reference sphere',fontsize=10)


def findIntercept(photonList,offset=0.0):
    numPhoton=len(photonList)
    stepPhoton=1
    intercept=[]
    for i in range(0,numPhoton,stepPhoton):
        dy=photonList[i].listY[-2]-photonList[i].listY[-3]
        dz=photonList[i].listZ[-2]-photonList[i].listZ[-3]
        intercept.append(-photonList[i].listY[-2]*dz/dy+photonList[i].listZ[-2])

    if numPhoton==1:
        intercept_mean=np.mean(np.array(intercept))
        print 'intercept', intercept_mean-offset, numPhoton
    else:
        intercept_mean=np.mean(np.array(intercept))
        intercept_std=np.std(np.array(intercept))
        print 'intercept', intercept_mean-offset, intercept_std, numPhoton
    return intercept_mean

def findClosestRay(photonList,verbose=1):
    numPhoton=len(photonList)
    stepPhoton=1
    minj=0
    mini=0
    minz=10.0
    for i in range(0,numPhoton,stepPhoton):
        #print photonList[i].xdir, photonList[i].ydir
        for j in range(len(photonList[i].listLayer)):
            #if photonList[i].listLayer[j]>100:
            #    print photonList[i].listLayer[j], photonList[i].listX[j], photonList[i].listY[j], photonList[i].listZ[j]
            if photonList[i].listLayer[j]==200:
                #print photonList[i].listX[j], photonList[i].listY[j], photonList[i].listZ[j]
                if np.fabs(photonList[i].listZ[j]-0.00290680501295)<minz:
                    minz=np.fabs(photonList[i].listZ[j]-0.00290680501295)
                    minj=j
                    mini=i
    if verbose>0:
        print 'closest ray:', photonList[mini].listX[minj], photonList[mini].listY[minj], photonList[mini].listZ[minj]
        print mini, minj
        print 'dz', minz
        for j in range(minj,len(photonList[mini].listLayer)):
            #print photonList[mini].listLayer[j], photonList[mini].listX[j], photonList[mini].listY[j], photonList[mini].listZ[j]
            print "%6g %12g %12g %12g %12g" % (photonList[mini].listLayer[j], photonList[mini].listX[j], photonList[mini].listY[j], photonList[mini].listZ[j], photonList[mini].opticZ[j])
        for j in range(minj+1,len(photonList[mini].listLayer)):
            print photonList[mini].listLayer[j], photonList[mini].opticZ[j]/np.sqrt((photonList[mini].listX[j]-photonList[mini].listX[j-1])**2.0+(photonList[mini].listY[j]-photonList[mini].listY[j-1])**2.0+(photonList[mini].listZ[j]-photonList[mini].listZ[j-1])**2.0)
    return mini,minj


def plotEvent(photonList,intercept=0.0,ii=-1):
    plotSlice=False
    numPhoton=len(photonList)
    startPhoton=0
    fig = plt.figure()
    ax=fig.add_subplot(111,aspect='equal')

    if ii==-1:
        plotPhoton= int(per * numPhoton)
        stepPhoton= int(1.0 / per)
        print plotPhoton, "out of", numPhoton
    else:
        startPhoton=ii
        numPhoton=ii+1
        stepPhoton=2

    #pdf = PdfPages(filename+'.pdf')

    plt.xlim([-1000,9000])
    plt.ylim([-5000,5000])

    print photonList[startPhoton].listLayer
    print photonList[startPhoton].listZ

    nn=0
    for i in range(startPhoton,numPhoton,stepPhoton):
       if plotSlice:
           if photonList[i].listX[8]<1000 and photonList[i].listX[8]>0:  #M1
               plt.plot(photonList[i].listZ[:-1],photonList[i].listY[:-1])
               nn+=1
       else:
           plt.plot(photonList[i].listZ[:-1],photonList[i].listY[:-1])
           print i, len(photonList[i].listZ[:-1])
           #print photonList[i].listLayer[-1], photonList[i].listLayer[-2], photonList[i].listLayer[-3], photonList[i].listLayer[-4]
    print nn

    plotOptics(opticsFile)
    plt.plot([-1000,9000],[0,0],':k')
    if intercept>0:
        plt.plot([intercept,intercept],[-5000,5000],':r')
        plt.text(intercept+40,4000,'exit pupil',fontsize=10)
        plt.plot(photonList[0].listZ[-3],photonList[0].listY[-3],'*m')
        plt.text(photonList[0].listZ[-3]+5,photonList[0].listY[-3],'image',fontsize=10)


    #pdf.savefig()
    plt.show()

def measureCrossSection(photonList):
    numPhoton=len(photonList)
    plotPhoton= int(per * numPhoton)
    stepPhoton= int(1.0 / per)
    print plotPhoton, "out of", numPhoton
    x=[]
    y=[]
    for i in range(0,numPhoton,stepPhoton):
        d=np.fabs(photonList[i].listZ[:-1])
        idx=d.argmin()
        if photonList[i].listZ[idx]<0 and photonList[i].listZ[idx-1]>0:
            fr=-photonList[i].listZ[idx]/(photonList[i].listZ[idx-1]-photonList[i].listZ[idx])
            dy=photonList[i].listY[idx]+(photonList[i].listY[idx-1]-photonList[i].listY[idx])*fr
            dx=photonList[i].listX[idx]+(photonList[i].listX[idx-1]-photonList[i].listX[idx])*fr
            x.append(dx*1e3) # in microns
            y.append(dy*1e3)
        else:
            print 'warning!', i, photonList[i].listZ[idx], photonList[i].listZ[idx-1]

    print 'std x, y', np.std(x), np.std(y)
    print 'min max x, y', np.min(x), np.max(x), np.min(y), np.max(y)

def plotIntercept(angles,intercepts):
    pdf = PdfPages('exitPupil.pdf')
    fig = plt.figure()
    ax=fig.add_subplot(111)
    plt.plot(angles,intercepts,'ko')
    ax.set_xlabel('off-axis angle (deg)')
    ax.set_ylabel('distance from AS (mm)')
    pdf.savefig()
    pdf.close()

def calculateOPL(X,Y,Z,n=[0]):
    opl=0.0
    if len(n)==1:
        ir=np.zeros(len(X))+1.0
    else:
        ir=n

    for i in range(len(X)-2):
        opl-=np.sqrt((X[i]-X[i+1])**2+(Y[i]-Y[i+1])**2+(Z[i]-Z[i+1])**2)*ir[i+1]
    i=len(X)-2
    opl-=np.sqrt((X[i]-X[i+1])**2+(Y[i]-Y[i+1])**2+(Z[i]-Z[i+1])**2)
    return opl

def checkOPDEvent(photonList):
    opl=[]
    distance=[]
    z=[]
    for i in range(len(photonList)):
        opl.append(calculateOPL(photonList[i].listX[7:21],photonList[i].listY[7:21],photonList[i].listZ[7:21])) #20km layer to Exit Pupil
        distance.append(np.sum(photonList[i].opticZ[8:20])-photonList[i].opticZ[20])
        z.append(photonList[i].listZ[7])
        if np.fabs(opl[i]/distance[i]-1)> 1e-7:
            print i, opl[i]/distance[i]

    l=0.0
    for j in range(8,21):
        l+=photonList[0].opticZ[j]*photonList[0].opticn[j]
        if j==20:
            l-=2*photonList[0].opticZ[j]*photonList[0].opticn[j]
        print '%d %11.4f %.3f %11.4f %11.4f' % (photonList[0].listLayer[j], photonList[0].opticZ[j], photonList[0].opticn[j], photonList[0].opticL[j], l)
    print 'bottom layer z', np.mean(z), np.std(z)

def verifyOPD(opdFile,photonList,ir=0,rays=100):
    #no n in opl
    opd=fits.getdata(opdFile)
    maxr=4180.0
    nn=len(opd[:,0])
    opdDiff=np.zeros((nn,nn))
    opdDiffEvent=np.zeros((nn,nn))
    opl=[]
    centralOPL=0.0
    checkR=1
    for i in range(len(photonList)):
        xx=int(photonList[i].listX[8]/maxr/2*(nn-1)+nn/2.0)
        yy=int(photonList[i].listY[8]/maxr/2*(nn-1)+nn/2.0)
        if rays==2:
            checkR=0
            if (xx==nn/2 and yy==nn/2) or (xx==127 and yy==nn-1):
                checkR=1
        #if checkR==0:
        #    opl.append(0)
        #    continue
        if ir==0:
            opl.append(calculateOPL(photonList[i].listX[7:21],photonList[i].listY[7:21],photonList[i].listZ[7:21])) #20km layer to Exit Pupil
        else:
            opl.append(calculateOPL(photonList[i].listX[7:21],photonList[i].listY[7:21],photonList[i].listZ[7:21],n=photonList[i].opticn[7:21])) #20km layer to Exit Pupil
        if (xx==nn/2 and yy==nn/2):
            centralOPL+=opl[i]
            print 'chief ray: %d %d %.6f %.6f %.6f\n' % (xx, yy, opl[i], photonList[i].opticL[20],opd[xx,yy])
        if (xx==127 and yy==254):
            print 'marginal ray: %d %d %.6f %.6f %.6f\n' % (xx, yy, opl[i], photonList[i].opticL[20],opd[xx,yy])
        if (xx==nn/2 and yy==10) or (xx==nn/2 and yy==40) or (xx==nn/2 and yy==nn/2) or (xx==127 and yy==nn-1):
            l=0.0
            print 'ray (%3d, %3d)' % (xx,yy)
            for j in range(8,21):
                l-=photonList[i].opticZ[j]*photonList[i].opticn[j]
                if j==20:
                    l-=2*photonList[i].opticZ[j]*photonList[i].opticn[j]
                print '%d %14.7f %14.7f %14.7f %14.7f %.8f %14.7f %14.7f' \
                % (photonList[i].listLayer[j], photonList[i].listX[j], photonList[i].listY[j], photonList[i].listZ[j], \
                        photonList[i].opticZ[j], photonList[i].opticn[j], photonList[i].opticL[j], l)
            print 'OPD %.5f\n' % (opd[xx,yy])




    print centralOPL
    #centralOPL=0.0
    #for i in range(len(photonList)):
    #    xx=int(photonList[i].listX[8]/maxr/2*(nn-1)+nn/2.0)
    #    yy=int(photonList[i].listY[8]/maxr/2*(nn-1)+nn/2.0)
    #    if opd[xx,yy]!=0:
    #        diffOPL=(opl[i]-centralOPL)-opd[xx,yy]
    #        opdDiff[xx,yy]=diffOPL/opd[xx,yy]
    #        #opdDiffEvent[xx,yy]=(photonList[i].opticL[20]-opd[xx,yy])/opd[xx,yy]
    #        #print xx, yy, opl[i]-centralOPL, opd[xx,yy]
    #fits.writeto('diff_'+opdFile, opdDiff, clobber=True)
    #fits.writeto('diffEvent_'+opdFile, opdDiffEvent, clobber=True)


def readZemaxRay(zemaxRay):
    zdata=np.zeros((13,6)) # x, y, z, pl, opl, cum opl
    fp=open(zemaxRay).readlines()
    for i in range(13):
        l=map(float, fp[i+1][13:].split())
        zdata[i,0]=l[0]
        zdata[i,1]=l[1]
        zdata[i,2]=l[2]
        zdata[i,3]=l[4]
        zdata[i,4]=l[5]
        zdata[i,5]=l[6]
    return zdata

def compareZemax(photonList,zemaxCR,zemaxMR):
    maxr=4180.0
    opl=[]
    centralOPL=0.0
    checkR=1
    nn=255
    zCR=readZemaxRay(zemaxCR)
    zMR=readZemaxRay(zemaxMR)
    pdf=PdfPages('position_diff.pdf')
    barLab=['M1','M2','M3','L1','L1E','L2','L2E','F','FE','L3','L3E','D','EP']
    ind=np.arange(13,0,-1)
    width=0.5
    plt.rcParams['font.size']=8
    for i in range(len(photonList)):
        xx=int(photonList[i].listX[8]/maxr/2*(nn-1)+nn/2.0)
        yy=int(photonList[i].listY[8]/maxr/2*(nn-1)+nn/2.0)

        l=0.0
        print 'ray (%3d, %3d)' % (xx,yy)
        print photonList[i].listX[8], photonList[i].listY[8], photonList[i].listZ[8]
        if xx==127 and yy==127:
            zd=zCR
            pt='Chief Ray'
            yl=[-0.5,0.5]
        else:
            zd=zMR
            pt='Marginal Ray'
            yl=[-4,4]
        k=0
        pd=np.zeros((13,6))
        for j in range(8,21):
            l-=photonList[i].opticZ[j]*photonList[i].opticn[j]
            print '%5s %d %14.7f %14.7f %14.7f %14.7f %14.7f %14.7f %.10f' \
            % (barLab[j-8], photonList[i].listLayer[j], photonList[i].listX[j]-zd[k,0], photonList[i].listY[j]-zd[k,1], -photonList[i].listZ[j]-zd[k,2], \
                    photonList[i].opticZ[j]-zd[k,3], photonList[i].opticn[j]*photonList[i].opticZ[j]-zd[k,4],\
                    photonList[i].opticL[j]+zd[k,5]+17.6094713729, photonList[i].opticn[j])
            pd[k,0]=photonList[i].listX[j]-zd[k,0]                          #dx
            pd[k,1]=photonList[i].listY[j]-zd[k,1]                          #dy
            pd[k,2]=-photonList[i].listZ[j]-zd[k,2]                         #dz
            pd[k,3]=photonList[i].opticZ[j]-zd[k,3]                         #pl
            pd[k,4]=photonList[i].opticn[j]*photonList[i].opticZ[j]-zd[k,4] #opl
            pd[k,5]=photonList[i].opticL[j]+zd[k,5]                         #cum_opl
            k+=1
        print ''
        #if xx==127 and yy==127:
        #    pd[:,1]+=1e-9

        pd=pd*1e6 # in nm
        for j in range(13):
            print '%5s %12.3e %12.3e %12.3e %12.3e' % (barLab[j],pd[j,4],pd[j,3],pd[j,2],pd[j,1])
        fig=plt.figure()
        fig.text(0.5,0.95,pt,ha='center',fontsize=10)
        ax=fig.add_subplot(141)
        ax.barh(ind, pd[:,4], width)
        ax.set_yticks(ind+width)
        ax.set_yticklabels(barLab)
        ax.set_xlabel(r'$\Delta$'+'OPL (nm)')
        #ax.set_xlim(yl)
        ax=fig.add_subplot(142)
        ax.barh(ind, pd[:,3], width)
        ax.set_yticks(ind+width)
        ax.set_yticklabels(barLab)
        ax.set_xlabel(r'$\Delta$'+'PL (nm)')
        #ax.set_xlim(yl)
        ax=fig.add_subplot(143)
        ax.barh(ind, pd[:,2], width)
        ax.set_yticks(ind+width)
        ax.set_yticklabels(barLab)
        ax.set_xlabel(r'$\Delta$'+'z (nm)')
        #ax.set_xlim(yl)
        ax=fig.add_subplot(144)
        ax.barh(ind, pd[:,1], width)
        ax.set_yticks(ind+width)
        ax.set_yticklabels(barLab)
        ax.set_xlabel(r'$\Delta$'+'y (nm)')
        #ax.set_xlim(yl)
        pdf.savefig()
    pdf.close()



eventFits = sys.argv[1]
per = float(sys.argv[2])
opticsFile = sys.argv[3]
print opticsFile
option = sys.argv[4]
#filename=eventFits.split('.fits')[0].split('output_')[1]
filename='opd'

#det_z=findHeight(opticsFile,'D')
if option=='intercept':
    photonList=readEvents(eventFits,verbose=0)
    offset=findHeight(opticsFile,'M1')
    intercept=findIntercept(photonList,offset)
elif option=='closest':
    photonList=readEvents(eventFits,verbose=0,addIntercept=1,addOP=1)
    mini,minj=findClosestRay(photonList)
elif option=='closestplot':
    photonList=readEvents(eventFits,verbose=0,addIntercept=1)
    mini,minj=findClosestRay(photonList,verbose=0)
    intercept=findIntercept([photonList[mini]],0.0)
    plotEvent(photonList,intercept,ii=mini)
elif option=='plotallintercept':
    offAxisAgnles=[0.001,0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    intercepts=[]
    for ang in offAxisAgnles:
        eventFits='output_ap5_'+str(ang)+'_d5_a0.fits.gz'
        photonList=readEvents(eventFits,verbose=0,addIntercept=1)
        mini,minj=findClosestRay(photonList,verbose=0)
        intercept=findIntercept([photonList[mini]],0.0)
        print ang, intercept
        intercepts.append(intercept)
    plotIntercept(offAxisAgnles,intercepts)
elif option=='checkopdevent':
    photonList=readEvents(eventFits,verbose=1,addOP=1)
    #checked: photon position diff = distance, less than 1e-6 errors
    #checkOPDEvent(photonList)
    opdFile = sys.argv[5]
    #verifyOPD(opdFile,photonList,ir=1,rays=2)
    #good to 1e-7 for 64-bit floating point number
elif option=='zemax':
    photonList=readEvents(eventFits,verbose=0,addOP=1)
    zemaxCR = sys.argv[5]
    zemaxMR = sys.argv[6]
    compareZemax(photonList,zemaxCR,zemaxMR)
elif option=='new':
    photonList=readEvents(eventFits,verbose=1,addOP=1)
    plotEvent(photonList)
elif option=='newopd':
    photonList=readEvents(eventFits,verbose=0)
    opdFile = sys.argv[5]
    verifyOPD(opdFile,photonList)
else:
    photonList=readEvents(eventFits,verbose=0,addIntercept=1)
    intercept=findIntercept(photonList,0.0)
    plotEvent(photonList,intercept)
#measureCrossSection(photonList)




