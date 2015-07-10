import subprocess, shutil
import math
import numpy as np
import optparse
import sys
sys.path.append("..")
from xyPosition import field2Sky, xyPositionRA, chipID
from optics import readOptics

def getSurfNum(surfName):
    num = []
    if surfName == 'M13':
        surfNames = ['M1','M3']
    #elif surfName in ['camera', 'Camera']:
    #    surfNames = ['L1','L1E','L2','L2E','F','FE','L3','L3E','D']
    else:
        surfNames = [surfName]
    for i, s in enumerate(surface):
        if s.name in surfNames:
            num.append(i)
    return num

def motionType(device,motion,d):
    if device == 'M2':
        num = [1]
    elif device == 'Camera':
        num = [3, 4, 5, 6, 7, 8, 9, 10, 11]
    else:
        print "Error"

    if motion == 'piston':
        typ = [[n, 5, -d] for n in num]
        fn = 'r1'
    elif motion == 'x-decenter':
        typ = [[n, 3, d] for n in num]
        fn = 'r2'
    elif motion == 'y-decenter':
        typ = [[n, 4, d] for n in num]
        fn = 'r3'
    elif motion == 'x-tilt':
        typ = axis2Euler('x',d, num)
        fn = 'r4'
    elif motion == 'y-tilt':
        typ = axis2Euler('y',d, num)
        fn = 'r5'
    else:
        print "Error"
    return typ, fn


def axis2Euler(axis,angle,surf):
    """
    R = | R11 R12 R13 |
        | R21 R22 R23 |
        | R31 R32 R33 |

    (x-convention)
    R11 = cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi)
    R12 = cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi)
    R13 = sin(psi)*sin(theta)
    R21 = -sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi)
    R22 = -sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi)
    R23 = cos(psi)*sin(theta)
    R31 = sin(theta)*sin(phi)
    R32 = -sin(theta)*cos(phi)
    R33 = cos(theta)


    x tilt
    Rx = | 1      0        0    |
         | 0   cos(rx)  sin(rx) |
         | 0  -sin(rx)  cos(rx) |

    y tilt
    Ry = | cos(ry)  0  -sin(ry) |
         |   0      1     0     |
         | sin(ry)  0   cos(ry) |
    """

    if axis == 'x':
        theta = -angle*math.pi/180
        phi = 0.0
        psi = 0.0
        cs = math.cos(angle*math.pi/180)
        ss = math.sin(angle*math.pi/180)
        R = np.array([[1,0,0],[0,cs,ss],[0,-ss,cs]])
    elif axis == 'y':
        theta = -angle*math.pi/180
        phi = math.pi/2
        psi = -math.pi/2
        cs = math.cos(angle*math.pi/180)
        ss = math.sin(angle*math.pi/180)
        R = np.array([[cs,0,-ss],[0,1,0],[ss,0,cs]])
    else:
        print 'Error'
    typ = [[surf[0], 0, phi], [surf[0], 1, psi], [surf[0], 2, theta]]
    p0 = np.array([[0], [0], [surface[surf[0]].height(0)]])
    for n in surf[1:]:
        p1 = np.array([[0], [0], [surface[n].height(0)]])
        p2 = np.dot(R,p1-p0) + p0
        dp = p2 - p1
        typ.append([n, 3, dp[0]])
        typ.append([n, 4, dp[1]])
        typ.append([n, 5, dp[2]])
        typ.append([n, 0, phi])
        typ.append([n, 1, psi])
        typ.append([n, 2, theta])
    return typ

def fieldPoint(i):
    degree=math.pi/180
    if i==1:
        fx, fy = 0.0, 0.0
    elif i==32:
        fx, fy = 1.185, 1.185
    elif i==33:
        fx, fy = -1.185, 1.185
    elif i==34:
        fx, fy = -1.185, -1.185
    elif i==35:
        fx, fy = 1.185, -1.185
    else:
        r = [0.379, 0.841, 1.237, 1.535, 1.708]
        theta = [0, 60, 120, 180, 240, 300]
        n = (i-2)/6
        m = (i-2)%6
        fx = r[n]*math.cos(theta[m]*degree)
        fy = r[n]*math.sin(theta[m]*degree)
    ra, dec = field2Sky(fx*degree,fy*degree)
    x, y = xyPositionRA(ra,dec)
    chip = chipID(x,y)
    return ra/degree, dec/degree, chip

def run(k,i=0,m=0,zernike=False,fea=None,surfName='M2',nollIdx=4,d=0.2):
    ra, dec, chip = fieldPoint(k)
    print k, chip, ra, dec
    inputPars = 'tmp'+str(i)+'_'+str(m)+'_'+str(nollIdx)+'.pars' if fea is None else 'tmp_'+fea[0:-4]+'.pars'
    comm = '../bin/raytrace < ' + inputPars
    if i == -1: #no perturbation
        fname = 'intrinsic_fld%d' % (k)
    else:
        if zernike:
            fname = '%s_z%d_%s_fld%d' % (surfName,nollIdx,d,k)
        elif fea is not None:
            fname = '%s_fld%d' % (fea[0:-4],k)
        else:
            device, motion = inputFile[i].split()[0:2]
            disp = inputFile[i+1].split()
            d = disp[m]
            typ, fn = motionType(device,motion,float(d))
            fname = '%s_%s_%s_fld%d' % (device,fn,d,k)

    pfile=open(inputPars,'w')
    pfile.write(open('raytrace_99999999_R22_S11_E000_opd0.pars').read())
    pfile.write('chipid %s\n' % (chip))
    pfile.write('opdfilename %s\n' % (fname))
    if i >= 0:
        if zernike:
            for n in getSurfNum(surfName):
                pfile.write('izernike %d %d %.12f\n' % (n, nollIdx-1, -d/1e3)) #d in microns
            if surfName == 'M13':
                pfile.write('izernikelink 2 0\n')
        elif fea is not None:
            for n in getSurfNum(surfName):
                pfile.write('fea %d %s\n' % (n, fea))
        else:
            for n, t, v in typ:
                pfile.write('body %d %d %.12f\n' % (n, t, v))

    pfile.write('object 0 %.10f %.10f 45 sed_0.50.txt 0 0 0 0 0 0 star none none\n' % (ra,dec))
    pfile.close()
    if subprocess.call(comm, shell=True) != 0:
        raise RuntimeError("Error running %s" % comm)



if __name__ == "__main__":
    inputFile = open('linearity_table_bending_short.txt').readlines()
    surface=readOptics('../data/lsst/optics_1.txt')
    parser = optparse.OptionParser()
    parser.add_option('-c',dest="col",default=0,type="int")
    parser.add_option('-r',dest="row",default=0,type="int")
    parser.add_option('-k',dest="field",default=1,type="int")
    parser.add_option('-s',dest="surfName",default='M2')
    parser.add_option('-n',dest="nollIdx",default=5,type="int")
    parser.add_option('-d',dest="d",default=0.1,type="float")
    parser.add_option('--all',action="store_true",dest="allfield",default=False)
    parser.add_option('--zernike',action="store_true",dest="zernike",default=False)
    parser.add_option('--fea',dest="fea",default=None)

    opt, remainder = parser.parse_args(sys.argv)
    if opt.allfield:
        for k in range(1,36):
            run(k,i=opt.row,m=opt.col,zernike=opt.zernike,surfName=opt.surfName,nollIdx=opt.nollIdx,d=opt.d,fea=opt.fea)
    else:
        run(opt.field,i=opt.row,m=opt.col,zernike=opt.zernike,surfName=opt.surfName,nollIdx=opt.nollIdx,d=opt.d,fea=opt.fea)
