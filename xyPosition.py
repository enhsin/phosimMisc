"""
   @file xyPosition.py
   @brief conversion between sky and image coordinates
   @brief Created by:
   @author En-Hsin Peng (Purdue)

Usage: python xyPosition.py ra 0 1.7   [convert sky (deg) to image (micron)]
       python xyPosition.py xy 0.0 306089.82   [convert image (micron) to sky (deg)]
       python xyPosition.py field 0 1.7   [convert field angle (deg) to image (micron)]


      x = cos(dec)sin(ra)
      y = sin(dec)
      z = cos(dec)cos(ra)
      p0 = (0, 0, 1)
      n = (0, 0, 1)
      l0 = (0, 0, 0)
      l = (x, y, z)
      d = (p0 -l0).n / (l.n)
        = 1/z
      p = d.l + l0
        = (x/z, y/z, 1)
      fx = atan2(1, x/z)
      fy = atan2(1, y/z)

"""

import sys
from math import sin, cos, pi, tan, atan2, atan
from scipy.optimize import fmin_tnc

DEGREE = pi/180.0
plateScale = 180000.0 # in data/lsst/central_wavelengths.txt
focalLength = plateScale/DEGREE

'''from trim:xyPosition'''
def xyPosition(alpha, delta, pointingRA, pointingDec):
    a = cos(delta)*cos(alpha - pointingRA)
    f = focalLength/(sin(pointingDec)*sin(delta) + a*cos(pointingDec))
    y = f*(cos(pointingDec)*sin(delta) - a*sin(pointingDec))
    x = f*cos(delta)*sin(alpha - pointingRA)
    return x, y

def xyPositionRA(alpha, delta, pointingRA = 0.0, pointingDec = 0.0, rotationAngle = 0.0):
    xp, yp = xyPosition(alpha, delta, pointingRA, pointingDec)
    x = (xp*cos(rotationAngle) + yp*sin(rotationAngle))
    y = (-xp*sin(rotationAngle) + yp*cos(rotationAngle))
    '''position in microns'''
    print x, y

def xyPositionField(angleX, angleY):
    x = focalLength*tan(angleX)
    y = focalLength*tan(angleY)
    print x, y

def field2Sky2(angleX, angleY):
    x = focalLength*tan(angleX)
    y = focalLength*tan(angleY)
    skyAngle(x, y)

def field2Sky(fx, fy):
    #x = cos(dec)sin(ra)
    #y = sin(dec)
    #z = cos(dec)cos(ra)
    if fx == 0:
        ra = fx
        dec = fy
    else:
        x = tan(fx)  #x/z
        y = tan(fy)  #y/z
        ra = fx
        dec = atan(y/(x/sin(ra)))
    print '%.10f %.10f' % (ra/DEGREE, dec/DEGREE)

def sky2Field(ra, dec):
    x = cos(dec)*sin(ra)
    y = sin(dec)
    z = cos(dec)*cos(ra)
    fx = atan2(x, z)
    fy = atan2(y, z)
    print '%.10f %.10f' % (fx/DEGREE, fy/DEGREE)

def solveXY(sky, ra0, dec0, x0, y0):
    xp, yp = xyPosition(sky[0], sky[1], ra0, dec0)
    return (xp-x0)**2+(yp-y0)**2

#Warning: use with caution
def skyAngle(x, y, rotationAngle = 0.0, pointingRA = 0.0, pointingDec = 0.0):
    x0 = (x*cos(rotationAngle) - y*sin(rotationAngle))
    y0 = (x*sin(rotationAngle) + y*cos(rotationAngle))
    ra0 = pointingRA
    dec0 = pointingDec
    b=[(-pi,pi),(-pi/2,pi/2)]
    mpar=fmin_tnc(solveXY, [0,pi/4], approx_grad=1, bounds=b, args=(ra0,dec0,x0,y0),messages=0, accuracy=1e-12)
    '''position in degrees'''
    print '%.10f %.10f' % (mpar[0][0]/DEGREE, mpar[0][1]/DEGREE)

if sys.argv[1] == 'ra':
    alpha=float(sys.argv[2])*DEGREE
    delta=float(sys.argv[3])*DEGREE
    xyPositionRA(alpha,delta)
elif sys.argv[1] == 'field':
    alpha=float(sys.argv[2])*DEGREE
    delta=float(sys.argv[3])*DEGREE
    xyPositionField(alpha,delta)
elif sys.argv[1] == 'fieldRA':
    alpha=float(sys.argv[2])*DEGREE
    delta=float(sys.argv[3])*DEGREE
    field2Sky(alpha,delta)
elif sys.argv[1] == 'fieldRA2':
    alpha=float(sys.argv[2])*DEGREE
    delta=float(sys.argv[3])*DEGREE
    field2Sky2(alpha,delta)
elif sys.argv[1] == 'raField':
    alpha=float(sys.argv[2])*DEGREE
    delta=float(sys.argv[3])*DEGREE
    sky2Field(alpha,delta)
elif sys.argv[1] == 'xy':
    x=float(sys.argv[2])
    y=float(sys.argv[3])
    skyAngle(x, y)

