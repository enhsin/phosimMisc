import sys
from astropy.io import fits
import numpy as np

def run(fName1,fName2,op):
    f1=fits.getdata(fName1)
    f2=fits.getdata(fName2)
    idx=np.nonzero(f1)
    if op == 'sub' or op == 'subimage':
        diff = f1 -f2
        diff2 = (f1 -f2)**2
        if op == 'subimage':
            fName = fName1[0:-8] + '_diff.fits.gz'
            fits.writeto(fName,diff,clobber=True)
        idx2=np.nonzero(diff)
    print 'mean: %14.6e %14.6e' % (np.mean(f1[idx]), np.mean(f2[idx]))
    print 'std:  %14.6e %14.6e' % (np.std(f1[idx]), np.std(f2[idx]))
    print 'max:  %14.6e %14.6e' % (np.max(f1[idx]), np.max(f2[idx]))
    print 'min:  %14.6e %14.6e' % (np.min(f1[idx]), np.min(f2[idx]))
    print 'mean diff: %11.3e' % np.mean(diff[idx])
    print 'std diff:  %11.3e' % np.std(diff[idx])
    print 'max diff:  %11.3e' % np.max(diff[idx])
    print 'min diff:  %11.3e' % np.min(diff[idx])
    print 'rms diff:  %11.3e' % np.mean(diff2[idx])**0.5
    print 'rms nonzero diff:  %11.3e %d' % (np.mean(diff2[idx2])**0.5, len(idx2[0]))


run(sys.argv[1],sys.argv[2],sys.argv[3])

