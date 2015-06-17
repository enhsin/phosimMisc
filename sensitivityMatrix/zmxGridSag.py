import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def readSag(fn):
    f = open(fn).readlines()
    nx, ny = map(int,f[0].split()[0:2])
    delx, dely = map(float,f[0].split()[2:4])
    z = np.zeros((nx,ny))
    z_x = np.zeros((nx,ny))
    z_y = np.zeros((nx,ny))
    z_xy = np.zeros((nx,ny))
    for j in range(ny):
        for i in range(nx):
            z[i,ny-j-1], z_x[i,ny-j-1], z_y[i,ny-j-1], z_xy[i,ny-j-1] = map(float,f[j*nx+i+1].split())
    print np.median(z), np.max(z), np.min(z)
    print np.median(z_x), np.max(z_x), np.min(z_x)
    print np.median(z_y), np.max(z_y), np.min(z_y)
    print np.median(z_xy), np.max(z_xy), np.min(z_xy)
    return z*1e6, z_x*1e6, z_y*1e6, z_xy*1e6

def setClim(imgplot,z,p=0.99):
    n=len(z)
    lidx=int(n*(1-p)/2)
    uidx=int(n*(1+p)/2-1)
    z.sort()
    minv=z[lidx]
    maxv=z[uidx]
    imgplot.set_clim([minv,maxv])

def plotSag(fn):
    """
    zemax manual p.336
    """
    z, z_x, z_y, z_xy = readSag(fn)

    idx=np.nonzero(z)
    plt.rcParams['font.size']=9
    fig=plt.figure()
    ax=fig.add_subplot(221)
    ax.set_axis_off()
    imgplot=plt.imshow(z)
    ax.set_title('z')
    setClim(imgplot,z[idx])
    cbar=plt.colorbar(imgplot,shrink=0.9,orientation='horizontal')
    for tick in cbar.ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(5)

    ax=fig.add_subplot(222)
    ax.set_axis_off()
    imgplot=plt.imshow(z_x)
    ax.set_title(r'$\partial z/\partial x$')
    setClim(imgplot,z_x[idx])
    cbar=plt.colorbar(imgplot,shrink=0.9,orientation='horizontal')
    for tick in cbar.ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(5)

    ax=fig.add_subplot(223)
    ax.set_axis_off()
    imgplot=plt.imshow(z_y)
    ax.set_title(r'$\partial z/\partial y$')
    setClim(imgplot,z_y[idx])
    cbar=plt.colorbar(imgplot,shrink=0.9,orientation='horizontal')
    for tick in cbar.ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(5)

    ax=fig.add_subplot(224)
    ax.set_axis_off()
    imgplot=plt.imshow(z_xy)
    ax.set_title(r'$\partial^{2}z/\partial x\partial y$')
    setClim(imgplot,z_xy[idx])
    cbar=plt.colorbar(imgplot,shrink=0.9,orientation='horizontal')
    for tick in cbar.ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(5)

    plt.savefig(fn+'.pdf')



plotSag(sys.argv[1])





