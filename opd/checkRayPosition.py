import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
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
        print self.listLayer[8], self.listX[8], self.listY[8], self.listZ[8], self.listZ[7], self.listY[19]

def readEvents(eventFits,verbose=1,addIntercept=0,addOP=0):
    plotEvery=False
    f= fits.open(eventFits)
    event=f[1].data
    numLine = len(event.field(3))
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
        ptb=r'$\Delta \rm{x} = \rm{x}_{\rm{phosim}}-\rm{x}_{\rm{zemax}}$'+', OPL = cumulative optical path length'
        if xx==127 and yy==127:
            zd=zCR
            pt='Difference in Ray Position, Chief Ray'
            yl=[-6,6]
        else:
            zd=zMR
            pt='Difference in Ray Position, General Ray'
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
        pleft=0.08
        pright=0.1
        pbottom=0.16
        pmargin=0.04
        ptop=0.1
        plw=(1-pright-pleft-3*pmargin)/4
        plh=(1-ptop-pbottom)
        fig=plt.figure()
        fig.text(0.5,0.95,pt,ha='center',fontsize=10)
        fig.text(0.5,0.05,ptb,ha='center',fontsize=10)
        ax=fig.add_axes([pleft+3*(plw+pmargin), pbottom, plw, plh])
        ax.barh(ind[:], pd[:,5], width)
        ax.set_yticks(ind+width)
        ax.set_yticklabels(barLab)
        ax.set_xlabel(r'$\Delta$'+'OPL (nm)')
        ax.set_xlim([-1,1])
        ax=fig.add_axes([pleft+2*(plw+pmargin), pbottom, plw, plh])
        ax.barh(ind, pd[:,2], width)
        ax.set_yticks(ind+width)
        ax.set_yticklabels(barLab)
        ax.set_xlabel(r'$\Delta$'+'z (nm)')
        ax.set_xlim(yl)
        ax=fig.add_axes([pleft+plw+pmargin, pbottom, plw, plh])
        ax.barh(ind, pd[:,1], width)
        ax.set_yticks(ind+width)
        ax.set_yticklabels(barLab)
        ax.set_xlabel(r'$\Delta$'+'y (nm)')
        ax.set_xlim(yl)
        ax=fig.add_axes([pleft, pbottom, plw, plh])
        ax.barh(ind, pd[:,0], width)
        ax.set_yticks(ind+width)
        ax.set_yticklabels(barLab)
        ax.set_xlabel(r'$\Delta$'+'x (nm)')
        ax.set_xlim(yl)
        pdf.savefig()
    pdf.close()

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
        conic=a[4]
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
        if innerRadius>0:
            v = np.linspace(0,innerRadius, 20)
            X = np.outer(v,np.cos(u))
            Y = np.outer(v,np.sin(u))
            R=np.sqrt(X**2 + Y**2)
            if radiusofcurv==0:
                Z=R-R+height
            else:
                Z=height-R**2/radiusofcurv/(1.0+np.sqrt(1.0-(conic+1.0)*R**2/radiusofcurv**2))
                for i in range(3,11):
                    Z-=a[i+2]*1e3*R**i
            plt.plot(Z,Y,'-',linewidth=1.5,color='#C0C0C0')

        if surfaceName[-1]!='E':
            if radiusofcurv==0:
                txtZ=height
            else:
                txtZ=height-outerRadius**2/radiusofcurv/(1.0+np.sqrt(1.0-(conic+1.0)*outerRadius**2/radiusofcurv**2))
            txtY=outerRadius+100
            if surfaceName in ['F','D']:
                txtY=-outerRadius-300
            plt.text(txtZ,txtY,surfaceName,fontsize=10)


    u = np.linspace(-np.pi/4, np.pi/4, 30)
    z = epR*np.cos(u) + cz
    y = epR*np.sin(u) + cy
    maxy = max(y)
    plt.plot(z, y,'--r')
    plt.text(min(z),maxy+100,'reference sphere',fontsize=10,ha='right')

    plt.plot([0,0],[-4180,4180],':k')
    plt.text(0,-4500,'entrance pupil',fontsize=10,ha='center')


def rotate(x,y,angle):
    x=np.array(x)
    y=np.array(y)
    angle=np.pi*angle/180
    xp=x*np.cos(angle)+y*np.sin(angle)
    yp=-x*np.sin(angle)+y*np.cos(angle)
    return xp, yp


def plotEvent(photonList,opticsFile):
    pdf = PdfPages('layout.pdf')
    fig = plt.figure()
    ax=fig.add_subplot(111,aspect='equal')
    ax.set_xlabel('z (mm)')
    ax.set_ylabel('x (mm)')

    plt.xlim([-1000,11000])
    plt.ylim([-5000,5000])
    plt.plot([-1000,11000],[0,0],'-k')
    plotOptics(opticsFile)
    clr=['b','g']
    txt=['CR','GR']
    zstart=10000
    for i in [0,1]:
        plt.plot(photonList[i].listZ[8:-1],photonList[i].listY[8:-1],color=clr[i])
        z=photonList[i].listZ[8]
        y=photonList[i].listY[8]
        dz=z-photonList[i].listZ[7]
        dy=y-photonList[i].listY[7]
        angleY=dy/np.sqrt(dz*dz+dy*dy)
        angleZ=dz/np.sqrt(dz*dz+dy*dy)
        d=y*angleY+(z-zstart)*angleZ
        y0=y-d*angleY
        z0=z-d*angleZ
        plt.plot([z,z0],[y,y0],color=clr[i])
        plt.text(z0+60,y0-50,txt[i],fontsize=10)
        y0=y-d*0.9*angleY
        z0=z-d*0.9*angleZ
        plt.plot(z0,y0,color=clr[i],marker='<',mec=clr[i])


    iz, iy = rotate([zstart,zstart],[-4000,2000],1.7)
    plt.plot(iz,iy,':k')
    opdy=-z*angleY/angleZ + y
    plt.plot([z,0],[y,opdy],':g')
    plt.text(-800,opdy,'opdx',fontsize=10)
    plt.plot(0,opdy,'g*',mec='g')
    plt.plot(photonList[0].listZ[-3],photonList[0].listY[-3],'y*',mec='y')
    pdf.savefig()
    pdf.close()


eventFits = sys.argv[1]
photonList=readEvents(eventFits,verbose=0,addOP=1)
cy=photonList[0].listY[-3]  #chief ray position on the image plane
cz=photonList[0].listZ[-3]
R0=2738.3087855563
epR=(R0**2+cy**2)**0.5

if len(sys.argv) == 4:
    zemaxCR = sys.argv[2]
    zemaxMR = sys.argv[3]
    compareZemax(photonList,zemaxCR,zemaxMR)
elif len(sys.argv) == 3:
    opticsFile = sys.argv[2]
    plotEvent(photonList,opticsFile)





