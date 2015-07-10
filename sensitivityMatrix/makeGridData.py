import numpy as np
import sys
sys.path.append("..")
from optics import readOptics

def zernike(n,d,x,y,rmax):
    z=0.
    z_x=0.
    z_y=0.
    z_xy=0.
    ds = d/1e3
    x2 = x**2; x3 = x**3; x4 = x**4; x5 = x**5
    y2 = y**2; y3 = y**3; y4 = y**4; y5 = y**5
    if n == 6:
        z = ds*np.sqrt(6.0)*(x*x-y*y)
    elif n == 7:
        z = ds*np.sqrt(8.0)*(3*x*x*y+3*y*y*y-2*y)
        z_x = ds*np.sqrt(8.0)*(6*x*y)/rmax
        z_y = ds*np.sqrt(8.0)*(3*x*x+9*y*y-2)/rmax
        z_xy = ds*np.sqrt(8.0)*(6*x)/rmax**2
    elif n == 19:
        z = ds*np.sqrt(12.0)*(15*x4*y+10*x2*y3-5*y5-12*x2*y+4*y3)
        z_x = ds*np.sqrt(12.0)*(60*x3*y+20*x*y3-24*x*y)/rmax
        z_y = ds*np.sqrt(12.0)*(15*x4+30*x2*y2-25*y4-12*x2+12*y2)/rmax
        z_xy = ds*np.sqrt(12.0)*(60*x3+60*x*y2-24*x)/rmax**2

    return z, z_x, z_y, z_xy



def makeGridData(noll,dx=None):
    n = 200
    surface = readOptics('../data/lsst/optics_1.txt')
    rmax = surface[0].rout
    if dx is None:
        dx = 2.0*rmax/(n-1)
        print dx
    maxy=0
    maxx=0

    pfile=open("M1_test_z"+str(noll)+"_new.DAT","w")
    pfile.write("%d %d %.8e %.8e\n" % (n,n,dx,dx))
    for j in range(n):
        y = - (j - n/2.0 + 0.5)*dx/rmax
        if y>maxy: maxy=y
        for i in range(n):
            x = (i - n/2.0 + 0.5)*dx/rmax
            if x>maxx: maxx=x
            if x*x+y*y > 1.0:
                pfile.write("%.8e %.8e %.8e %.8e\n"  % (0,0,0,0))
            else:
                z, z_x, z_y, z_xy = zernike(noll,0.2,x,y,rmax)
                pfile.write("%.8e %.8e %.8e %.8e\n"  % (z,z_x,z_y,z_xy))

    pfile.close()
    print '%.3e %.3e' % (maxy, maxx)

def makeFeaData(noll):
    n = 150
    surface = readOptics('../data/lsst/optics_1.txt')
    rmax = surface[0].rout
    rmin = surface[0].rin
    dx = 2.0*rmax/(n-1)
    pfile=open("M1_test_z"+str(noll)+".dat","w")
    for j in range(n):
        y = - (j - n/2.0 + 0.5)*dx
        yi = y/rmax
        for i in range(n):
            x = (i - n/2.0 + 0.5)*dx
            xi = x/rmax
            r = np.sqrt(x*x+y*y)
            if r <= rmax and r>= rmin:
                z, z_x, z_y, z_xy = zernike(noll,0.2,xi,yi,rmax)
                pfile.write("%15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e\n"
                           % (x,y,0,0,0,-z,0,0,0))
    pfile.close()


#makeGridData(19,dx=42.0)
#makeGridData(19)
makeFeaData(19)
