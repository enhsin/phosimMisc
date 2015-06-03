import subprocess, shutil
import math
import sys
sys.path.append("..")
from xyPosition import field2Sky, xyPositionRA, chipID

def deviceNum(device):
    if device == 'M2':
        num = [1]
    elif device == 'camera':
        num = [3, 4, 5, 6, 7, 8, 9, 10]
    return num

def motionType(motion,d):
    if motion == 'piston':
        typ = [[5, -d]]
        fn = 'r1'
    elif motion == 'x-decenter':
        typ = [[3, d]]
        fn = 'r2'
    elif motion == 'y-decenter':
        typ = [[4, d]]
        fn = 'r3'
    elif motion == 'x-tilt':  # z = x tan(theta_x) + y tan(theta_y)
        phi, theta, psi = axis2Euler('x',d)
        typ = [[0, phi], [1, psi], [2, theta]]
        fn = 'r4'
    elif motion == 'y-tilt':
        phi, theta, psi = axis2Euler('y',d)
        typ = [[0, phi], [1, psi], [2, theta]]
        fn = 'r5'
    else:
        print "Error"
    return typ, fn


def axis2Euler(axis,angle):
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
    elif axis == 'y':
        theta = -angle*math.pi/180
        phi = math.pi/2
        psi = -math.pi/2
    else:
        print 'Error'
    return phi, theta, psi


def fieldPoint(i):
    if i==1:
        fx, fy = 0.0, 0.0
    else:
        r = [0.379, 0.841, 1.237, 1.535, 1.708]
        theta = [0, 60, 120, 180, 240, 300]
        n = (i-2)/6
        m = (i-2)%6
        fx = r[n]*math.cos(theta[m]*math.pi/180)*math.pi/180
        fy = r[n]*math.sin(theta[m]*math.pi/180)*math.pi/180
    ra, dec = field2Sky(fx,fy)
    x, y = xyPositionRA(ra,dec)
    chip = chipID(x,y)
    return ra, dec, chip

def run(k,i):
    ra, dec, chip = fieldPoint(k)
    print k, chip, ra, dec
    comm = '../bin/raytrace < tmp.pars'
    device, motion = inputFile[i].split()[0:2]
    disp = inputFile[i+1].split()
    for d in [disp[0], disp[4]]:
        typ, fn = motionType(motion,float(d))
        opdName = 'opd_%s_%s_%s_fld%d' % (device,fn,d,k)
        pfile=open('tmp.pars','w')
        pfile.write(open('raytrace_99999999_R22_S11_E000_opd0.pars').read())
        pfile.write('chipid %s\n' % (chip))
        for n in deviceNum(device):
            for t, v in typ:
                pfile.write('body %d %d %f\n' % (n, t, v))
        pfile.write('object 0 %.10f %.10f 45 sed_0.50.txt 0 0 0 0 0 0 star none none\n' % (ra,dec))
        pfile.close()
        if subprocess.call(comm, shell=True) != 0:
            raise RuntimeError("Error running %s" % comm)
        else:
            shutil.move('opd.fits.gz',opdName+'.fits.gz')
            shutil.move('opd_20415.txt',opdName+'_20415.txt')
            shutil.move('opd_0.txt',opdName+'_0.txt')
            shutil.move('lsst_e_99999999_'+chip+'_E000.fits','lsst_e_99999999_'+chip+'_E000_fld'+str(k)+'.fits')


inputFile = open('linearity_table_bending_short.txt').readlines()
for k in range(1,32):
    ra, dec, chip = fieldPoint(k)
    print k, chip, ra, dec

#for k in range(1,32):
#    run(k,0)
#    run(k,2)
#    run(k,4)
#run(int(sys.argv[1]),int(sys.argv[2]))
