import sys
from astropy.io import fits
import numpy as np

def run(fName1,fName2,op):
    f1=fits.getdata(fName1)
    f2=fits.getdata(fName2)
    idx=np.nonzero(f1)
    if op == 'sub':
        diff = f1 -f2
    print 'mean: %14.6e %14.6e' % (np.mean(f1[idx]), np.mean(f2[idx]))
    print 'std:  %14.6e %14.6e' % (np.std(f1[idx]), np.std(f2[idx]))
    print 'max:  %14.6e %14.6e' % (np.max(f1[idx]), np.max(f2[idx]))
    print 'min:  %14.6e %14.6e' % (np.min(f1[idx]), np.min(f2[idx]))
    print 'mean diff: %11.3e' % np.mean(diff[idx])
    print 'std diff:  %11.3e' % np.std(diff[idx])
    print 'max diff:  %11.3e' % np.max(diff[idx])
    print 'min diff:  %11.3e' % np.min(diff[idx])


run(sys.argv[1],sys.argv[2],sys.argv[3])

