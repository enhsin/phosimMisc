from astropy.io import fits
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plotZMX(pfn,zfn):
    tx='field angle = 1.7 deg'
    zemax0=np.loadtxt(zfn)
    zemax=zemax0[:-1,1:]
    idx=np.nonzero(zemax==0)
    idx2=np.nonzero(zemax)
    print sum(zemax0[:,0]), sum(zemax0[-1,:])
    n=len(idx2[0])
    p=0.98
    lidx=int(n*(1-p)/2)
    uidx=int(n*(1+p)/2-1)
    simage=zemax[idx2]
    simage.sort()
    minv=simage[lidx]
    maxv=simage[uidx]

    plt.rcParams['font.size']=9
    fig=plt.figure()
    fig.text(0.5,0.95,tx,ha='center')
    ax=fig.add_subplot(221)
    ax.set_axis_off()
    image2=np.zeros((255,255))+minv-2*(maxv-minv)
    image2[idx2]=zemax[idx2]
    imgplot=plt.imshow(image2)
    imgplot.set_clim([minv,maxv])
    ax.set_title('zemax')
    cbar=plt.colorbar(imgplot,shrink=0.9,orientation='horizontal')
    cbar.set_label('#')
    for tick in cbar.ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(5)

    phosim1=fits.getdata(pfn)*1000/0.77
    phosim=np.transpose(phosim1)
    phosim=phosim[::-1,:]
    ax=fig.add_subplot(222)
    ax.set_axis_off()
    image2=np.zeros((255,255))+minv-2*(maxv-minv)
    image2[idx2]=phosim[idx2]
    imgplot=plt.imshow(image2)  #in mm/770 nm
    imgplot.set_clim([minv,maxv])
    ax.set_title('phosim')
    cbar=plt.colorbar(imgplot,shrink=0.9,orientation='horizontal')
    cbar.set_label('#')
    for tick in cbar.ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(5)

    ax=fig.add_subplot(223)
    ax.set_axis_off()
    diff=(phosim-zemax)
    fits.writeto('diff.fits.gz',diff,clobber=True)
    simage=diff[idx2]
    simage.sort()
    minv=simage[lidx]
    maxv=simage[uidx]
    image2=np.zeros((255,255))+minv-2*(maxv-minv)
    image2[idx2]=diff[idx2]
    imgplot=plt.imshow(image2*770)
    imgplot.set_clim([minv*770,maxv*770])
    print np.median(np.fabs(diff[idx2])), np.min(diff[idx2]), np.max(diff[idx2])
    ax.set_title('phosim-zemax')
    cbar=plt.colorbar(imgplot,shrink=0.9,orientation='horizontal')
    cbar.set_label('nm')
    for tick in cbar.ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(5)

    ax=fig.add_subplot(224)
    ax.set_axis_off()
    zemax[idx]=1.0
    diff=diff/zemax
    simage=diff[idx2]
    simage.sort()
    minv=simage[lidx]
    maxv=simage[uidx]
    image2=np.zeros((255,255))+minv-2*(maxv-minv)
    image2[idx2]=diff[idx2]
    imgplot=plt.imshow(image2*100)
    imgplot.set_clim([minv*100,maxv*100])

    print np.mean(np.fabs(diff[idx2])), np.min(diff[idx2]), np.max(diff[idx2])
    ax.set_title('(phosim-zemax)/zemax')
    cbar=plt.colorbar(imgplot,shrink=0.9,orientation='horizontal')
    cbar.set_label('%')
    for tick in cbar.ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(5)


    plt.savefig('opd.pdf')



plotZMX(sys.argv[1],sys.argv[2])
