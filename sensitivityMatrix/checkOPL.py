import sys
import numpy as np

def run(fName1,fName2,col):
    f1=np.loadtxt(fName1)
    f2=np.loadtxt(fName2)
    for i in range(len(f1[:,0])):
        print '%17.10f %17.10f %17.2f' % (f1[i,col],f2[i,col],(f1[i,col]-f2[i,col])*1e6)

def compareZemax(fzemax,fphosim,col):
    f2=np.loadtxt(fphosim)
    f1=np.zeros((13,7))
    fp=open(fzemax).readlines()
    for i in range(13):
        l=map(float, fp[i+1][13:].split())
        f1[i,1]=l[0]  #x
        f1[i,2]=l[1]  #y
        f1[i,3]=-l[2]  #z
        f1[i,4]=l[4]  #pl
        f1[i,5]=l[3]  #n
        f1[i,6]=-l[6]  #cum opl
    for i in range(len(f1[:,0])):
        print '%17.10f %17.10f %17.2f' % (f1[i,col],f2[i,col],(f1[i,col]-f2[i,col])*1e6)


if 'zmx' in sys.argv:
    compareZemax(sys.argv[1],sys.argv[2],int(sys.argv[3]))
else:
    run(sys.argv[1],sys.argv[2],int(sys.argv[3]))


