import sys
from astropy.io import fits
import numpy as np

def run(fName1,fName2,op):
    f1=fits.getdata(fName1)
    f2=fits.getdata(fName2)
    idx=np.nonzero(f1)
    if op == 'sub':
        diff = f1 -f2
    print 'mean:', np.mean(f1[idx])
    print 'std:', np.std(f1[idx])
    print 'mean diff:', np.mean(diff[idx])
    print 'std diff:', np.std(diff[idx])


run(sys.argv[1],sys.argv[2],sys.argv[3])

