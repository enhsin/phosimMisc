import numpy as np
import os, glob
from astropy.io import fits

def coadd(files, outfile):
    hdr = fits.getheader(files[0])
    images = [fits.getdata(f) for f in files]
    data = np.median(images, axis=0)
    fits.writeto(outfile, data, hdr, clobber=True)
    os.system('gzip -f '+outfile)

def imageList(root, cal, imgType, fid):
    obsids=[calObs[cal]+str(s) for s in range(10)]
    if imgType == 'e':
        f=[root+cal+'/eimage/'+obs+'/lsst_e_'+obs+'_'+fid+'.fits.gz'  for obs in obsids]
    else:
        f=[root+cal+'/raw/'+obs+'/lsst_a_'+obs+'_'+fid+'.fits.gz'  for obs in obsids]
    return f


def coaddAll(root, cal, imgType):
    if imgType == 'e':
        fileIDs = ['f'+f.split('/')[-1].split('.')[0].split('f')[-1] \
                for f in glob.glob(root+cal+'/eimage/'+calObs[cal]+'0/*fits.gz')]
        for fid in fileIDs:
            f = imageList(root, cal, 'e', fid)
            #print cal+'/eimage/lsst_e_'+calObs[cal][0:-1]+'10_'+fid+'.fits'
            coadd(f, cal+'/eimage/lsst_e_'+calObs[cal][0:-1]+'10_'+fid+'.fits')
    else:
        fileIDs = ['f'+f.split('/')[-1].split('.')[0].split('f')[-1] \
                for f in glob.glob(root+cal+'/raw/'+calObs[cal]+'0/*fits.gz')]
        for fid in fileIDs:
            f = imageList(root, cal, 'a', fid)
            #print cal+'/raw/lsst_a_'+calObs[cal][0:-1]+'10_'+fid+'.fits'
            coadd(f, cal+'/raw/lsst_a_'+calObs[cal][0:-1]+'10_'+fid+'.fits')



calObs={'bias': '9999000', 'dark': '9999100',
        'flat_u': '9999200', 'flat_g': '9999300',
        'flat_r': '9999400', 'flat_i': '9999500',
            'flat_z': '9999600', 'flat_y': '9999700'}
root='/depot/lsst/phosim/data/'
coaddAll(root, 'dark', 'e')
coaddAll(root, 'dark', 'a')
coaddAll(root, 'bias', 'e')
coaddAll(root, 'bias', 'a')
