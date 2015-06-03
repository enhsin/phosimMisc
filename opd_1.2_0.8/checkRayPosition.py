import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages


def readZemaxRay(zemaxRay):
    zdata=np.zeros((13,6))
    fp=open(zemaxRay).readlines()
    for i in range(13):
        l=map(float, fp[i+1][13:].split())
        zdata[i,0]=l[0]  #x
        zdata[i,1]=l[1]  #y
        zdata[i,2]=l[2]  #z
        zdata[i,3]=l[3]  #n
        zdata[i,4]=l[4]  #pl
        zdata[i,5]=l[6]  #cum opl
    return zdata

def readPhosimRay(phosimRay):
    zdata=np.zeros((13,6))
    data=np.loadtxt(phosimRay)
    zdata[:,0]=data[:,1]  #x
    zdata[:,1]=data[:,2]  #y
    zdata[:,2]=-data[:,3] #z
    zdata[:,3]=data[:,5]  #n
    zdata[:,4]=data[:,4]  #pl
    zdata[:,5]=-data[:,6] #cum opl
    return zdata

def compareZemax(phosimR,zemaxR):
    zR=readZemaxRay(zemaxR)
    pR=readPhosimRay(phosimR)
    pdf=PdfPages('position.pdf')
    barLab=['M1','M2','M3','L1','L1E','L2','L2E','F','FE','L3','L3E','D','EP']
    ind=np.arange(13,0,-1)
    width=0.5
    plt.rcParams['font.size']=8
    pt='Difference in Ray Position'
    ptb=r'$\Delta \rm{x} = \rm{x}_{\rm{phosim}}-\rm{x}_{\rm{zemax}}$'+', OPL = cumulative optical path length'
    yl=[-6,6]
    pd=(pR-zR)*1e6 # in nm
    print '%4s  %7s %7s %7s %8s %7s %7s' % ('surf','dx(nm)','dy(nm)','dz(nm)','dn(1e-6)','dpl(nm)','dcumopl')
    for i in range(13):
        print '%4s  %7.2f %7.2f %7.2f %8.2f %7.2f %7.2f' % (barLab[i], pd[i,0], pd[i,1], pd[i,2], pd[i,3], pd[i,4], pd[i,5])
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
    #ax.set_xlim([-1,1])
    ax=fig.add_axes([pleft+2*(plw+pmargin), pbottom, plw, plh])
    ax.barh(ind, pd[:,2], width)
    ax.set_yticks(ind+width)
    ax.set_yticklabels(barLab)
    ax.set_xlabel(r'$\Delta$'+'z (nm)')
    #ax.set_xlim(yl)
    ax=fig.add_axes([pleft+plw+pmargin, pbottom, plw, plh])
    ax.barh(ind, pd[:,1], width)
    ax.set_yticks(ind+width)
    ax.set_yticklabels(barLab)
    ax.set_xlabel(r'$\Delta$'+'y (nm)')
    #ax.set_xlim(yl)
    ax=fig.add_axes([pleft, pbottom, plw, plh])
    ax.barh(ind, pd[:,0], width)
    ax.set_yticks(ind+width)
    ax.set_yticklabels(barLab)
    ax.set_xlabel(r'$\Delta$'+'x (nm)')
    #ax.set_xlim(yl)
    pdf.savefig()
    pdf.close()


phosimR = sys.argv[1]
zemaxR = sys.argv[2]
compareZemax(phosimR,zemaxR)





