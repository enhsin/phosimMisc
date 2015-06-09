from astropy.io import fits
import numpy as np
import os

def rebin(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)

def resize(img,xNew,yNew):
    xOrg=np.arange(len(img[:,0]))
    yOrg=np.arange(len(img[0,:]))
    #f=interpolate.interp2d(xOrg, yOrg, img, kind='linear', copy=False)
    #f=interpolate.interp2d(xOrg, yOrg, img, kind='linear')
    f=interpolate.RectBivariateSpline(xOrg, yOrg, img, kx=1, ky=1)
    imgNew=f(xNew,yNew)
    return imgNew

def mosaicEimage(obsids):
    instrDir='/scratch/hansen/e/epeng/sims_phosim/data/lsst'
    #instrDir='/home/abby/PROGRAM/phosim/data/lsst'
    fp=open(instrDir+"/focalplanelayout.txt").readlines()
    chipID=[]
    xpos=[]
    ypos=[]
    numpx=[]
    numpy=[]
    x1=[]
    x2=[]
    y1=[]
    y2=[]
    pixsize=10.0
    sq=640*10
    psize=100.0
    for line in fp:
        lstr=line.split()
        if "Group0" in line:
            chipID.append(lstr[0])
            xpos.append(float(lstr[1]))
            ypos.append(float(lstr[2]))
            numpx.append(float(lstr[4]))
            numpy.append(float(lstr[5]))
            x1.append(0)
            x2.append(float(lstr[4])-1)
            y1.append(0)
            y2.append(float(lstr[5])-1)

    for obsid in obsids:
        fld='/apps/group/lsst/cal_v3.4/telrot/'+obsid
        #fld='.'
        fullimage=np.zeros((sq,sq))
        for i in range(len(chipID)):
            #if 'R22' not in chipID[i]:
            #    continue
            f=fits.open(fld+'/lsst_e_'+obsid+'_f0_'+chipID[i]+'_E000.fits.gz')
            img=f[0].data
            ny=len(img[0,:])
            nx=len(img[:,0])
            x0=int(np.round(xpos[i]/psize+x1[i]*pixsize/psize-numpx[i]*pixsize/psize/2.0+sq/2))
            y0=int(np.round(ypos[i]/psize+y1[i]*pixsize/psize-numpy[i]*pixsize/psize/2.0+sq/2))
            imgsmall=np.transpose(rebin(img[1:nx-1,:],[407,400]))
            xc=[x0,x0+400]
            yc=[sq-y0-407,sq-y0]
            x0c=[0,400]
            y0c=[0,407]
            if xc[0]<0:
                xc[0]=0
                x0c[0]=-x0
            if xc[1]>sq:
                xc[1]=sq
                x0c[1]=sq-x0
            if yc[0]<0:
                yc[0]=0
                y0c[1]=sq-y0
            if y0c[1]>sq:
                yc[1]=sq
                y0c[0]=-y0
            fullimage[xc[0]:xc[1],yc[0]:yc[1]]+=imgsmall[x0c[0]:x0c[1],y0c[0]:y0c[1]][:,::-1]

        fits.writeto('eimage_'+obsid+'.fits', fullimage, clobber=True)
        os.system('gzip -f eimage_'+obsid+'.fits')

def plotEimage(obsids):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    pdf=PdfPages('flat_u_rot.pdf')
    for obsid in obsids:
        angle=(float(obsid)-99992000)/10
        fullimage=fits.getdata('eimage_'+obsid+'.fits.gz')
        sq=len(fullimage[:,0])
        fig=plt.figure()
        ax=fig.add_subplot(111)
        ax.set_axis_off()
        smallimage=rebin(fullimage,[sq/5,sq/5])
        #maxlum=np.max(smallimage)*0.95
        #minlum=maxlum*0.1
        #idx=np.nonzero(smallimage>maxlum)
        #smallimage[idx]=maxlum
        #idx=np.nonzero(smallimage<minlum)
        #smallimage[idx]=minlum
        #scaledImage=(smallimage-minlum)/(maxlum-minlum)
        #imgplot=plt.imshow(scaledImage)
        if obsid=='99992100':
            smallimage0=smallimage
            fullimage0=fullimage
            imgplot=plt.imshow(smallimage)
            bad=np.nonzero(smallimage0==0)
            smallimage0[bad]=1.0
            bad=np.nonzero(fullimage0==0)
            fullimage0[bad]=1.0
            ax.set_title('%.0f deg' % angle)
        else:
            imgplot=plt.imshow(smallimage/smallimage0)
            fits.writeto('eimage_'+obsid+'_div.fits', fullimage/fullimage0, clobber=True)
            os.system('gzip -f eimage_'+obsid+'_div.fits')
            imgplot.set_clim([0.988, 1.012])
            ax.set_title('%.0f deg / 10 deg' % angle)
        cbar=plt.colorbar(imgplot,shrink=0.9,orientation='vertical')
        #imgplot.set_cmap('gist_gray')

        pdf.savefig()
    pdf.close()


#mosaicEimage(['99992300','99992400'])
#mosaicEimage(['99992100','99992200','99992300','99992400','99992500'])
plotEimage(['99992100','99992200','99992300','99992400','99992500'])
