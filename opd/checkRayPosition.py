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
        if verbose>0:
            if i%numLine20==0:
                print j,'/20'
                j+=1
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
    pdf=PdfPages('position.pdf')
    barLab=['M1','M2','M3','L1','L1E','L2','L2E','F','FE','L3','L3E','D','EP']
    ind=np.arange(13,0,-1)
    width=0.5
    plt.rcParams['font.size']=8
    for i in range(len(photonList)):
        xx=int(photonList[i].listX[8]/maxr/2*(nn-1)+nn/2.0)
        yy=int(photonList[i].listY[8]/maxr/2*(nn-1)+nn/2.0)

        l=0.0
        print 'ray (%3d, %3d)' % (xx,yy)
        if xx==127 and yy==127:
            zd=zCR
            pt='Difference between PhoSim and Zemax, Chief Ray'
            yl=[-6,6]
        else:
            zd=zMR
            #pt='Marginal Ray'
            pt='Difference between PhoSim and Zemax, General Ray'
            yl=[-6,6]
        k=0
        pd=np.zeros((13,6))
        print '%9s %14s %14s %14s %14s %14s %12s' % ('Surface','dx(mm)','dy(mm)','dz(mm)','dd(mm)','dopl(mm)','n') 
        for j in range(8,21):
            l-=photonList[i].opticZ[j]*photonList[i].opticn[j]
            print '%9s %14.7f %14.7f %14.7f %14.7f %14.7f %.10f' \
            % (barLab[j-8], photonList[i].listX[j]-zd[k,0], photonList[i].listY[j]-zd[k,1], -photonList[i].listZ[j]-zd[k,2], \
                    photonList[i].opticZ[j]-zd[k,3], \
                    photonList[i].opticL[j]+zd[k,5], photonList[i].opticn[j])
            pd[k,0]=photonList[i].listX[j]-zd[k,0]                          #dx
            pd[k,1]=photonList[i].listY[j]-zd[k,1]                          #dy
            pd[k,2]=-photonList[i].listZ[j]-zd[k,2]                         #dz
            pd[k,3]=photonList[i].opticZ[j]-zd[k,3]                         #pl
            pd[k,4]=photonList[i].opticn[j]*photonList[i].opticZ[j]-zd[k,4] #opl
            pd[k,5]=photonList[i].opticL[j]+zd[k,5]                         #cum_opl
            k+=1
        print ''

        pd=pd*1e6 # in nm
        fig=plt.figure()
        fig.text(0.5,0.95,pt,ha='center',fontsize=10)
        ax=fig.add_subplot(144)
        ax.barh(ind[:], pd[:,5], width)
        ax.set_yticks(ind+width)
        ax.set_yticklabels(barLab)
        ax.set_xlabel('OPL (nm)')
        ax.set_xlim([-1,1])
        ax=fig.add_subplot(143)
        ax.barh(ind, pd[:,2], width)
        ax.set_yticks(ind+width)
        ax.set_yticklabels(barLab)
        ax.set_xlabel(r'$\Delta$'+'z (nm)')
        ax.set_xlim(yl)
        ax=fig.add_subplot(142)
        ax.barh(ind, pd[:,1], width)
        ax.set_yticks(ind+width)
        ax.set_yticklabels(barLab)
        ax.set_xlabel(r'$\Delta$'+'y (nm)')
        ax.set_xlim(yl)
        ax=fig.add_subplot(141)
        ax.barh(ind, pd[:,0], width)
        ax.set_yticks(ind+width)
        ax.set_yticklabels(barLab)
        ax.set_xlabel(r'$\Delta$'+'x (nm)')
        ax.set_xlim(yl)
        pdf.savefig()
    pdf.close()



eventFits = sys.argv[1]
zemaxCR = sys.argv[2]
zemaxMR = sys.argv[3]

photonList=readEvents(eventFits,verbose=0,addOP=1)
compareZemax(photonList,zemaxCR,zemaxMR)




