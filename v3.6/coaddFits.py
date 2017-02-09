import numpy as np
import os, glob, sys
from astropy.io import fits

def coadd(files, outputDir, log=None):
    outfile = outputDir+files[0][0:13]+'20'+files[0][15:-9]+'0.fits'
    hdr = fits.getheader(files[0])
    images = [fits.getdata(f) for f in files]
    data = np.median(images, axis=0)
    fits.writeto(outfile, data, hdr, clobber=True)
    os.system('gzip -f '+outfile)
    if log is not None:
        sum = np.array([im.sum() for im in images])
        m = sum.mean()
        s = sum.std()
        r = s/m
        log.write('%s %.0f %.0f %.5f\n' % (files[0][16:-13],m,s,r))
        if r > 0.01:
            print '%s %.0f %.0f %.5f' % (files[0][16:-13],m,s,r)

def coaddAll(intputDir, outputDir, cal):
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    if not os.path.exists(outputDir+'/raw'):
        os.makedirs(outputDir+'/raw')
    if not os.path.exists(outputDir+'/eimage'):
        os.makedirs(outputDir+'/eimage')

    cur=os.path.split(os.path.abspath(__file__))[0]
    logfile=open(cal+'.dat','w')
    os.chdir(inputDir)
    for f in glob.glob('*'+calObs[cal]+'0*R*_S*E000.tar'):
    #for f in glob.glob('*'+calObs[cal]+'0*R22_S11*E000.tar'):
        fid=f[14:25]
        for i in range(10):
            os.system('tar xf lsst_'+calObs[cal]+str(i)+'_'+fid+'E000.tar')
            os.system('tar xf lsst_'+calObs[cal]+str(i)+'_'+fid+'E001.tar')

        images = glob.glob('lsst_e_'+calObs[cal]+'*.fits.gz')
        if len(images)==20:
            coadd(images, outputDir+'/eimage/')
        else:
            print images
            print 'eimage errors', f
        for im in images:
            os.remove(im)

        for aid in ['C00','C01','C02','C03','C04','C05','C06','C07','C10','C11','C12','C13','C14','C15','C16','C17']:
            images = glob.glob('lsst_a_'+calObs[cal]+'*'+aid+'*.fits.gz')
            if len(images)==20:
                coadd(images, outputDir+'/raw/',logfile)
            else:
                print 'raw image errors', f, aid
            for im in images:
                os.remove(im)
    os.chdir(cur)
    logfile.close()


calObs={'bias': '9999000', 'dark': '9999100',
        'flat_u': '9999200', 'flat_g': '9999300',
        'flat_r': '9999400', 'flat_i': '9999500',
            'flat_z': '9999600', 'flat_y': '9999700'}
inputDir=os.path.split(os.path.abspath(__file__))[0]+'/output'
print inputDir
sys.exit()
outputDir='/depot/lsst/phosim/data/v3.6/'
#coaddAll(inputDir, outputDir+'dark', 'dark')
#coaddAll(inputDir, outputDir+'bias', 'bias')
for f in ['u','g','r','i','z','y']:
    coaddAll(inputDir, outputDir+'flat_'+f, 'flat_'+f)
