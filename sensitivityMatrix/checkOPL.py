import sys
import numpy as np

def run(fName1,fName2,col):
    f1=np.loadtxt(fName1)
    f2=np.loadtxt(fName2)
    for i in range(len(f1[:,0])):
        print '%17.10f %17.10f %17.10f' % (f1[i,col],f2[i,col],f1[i,col]-f2[i,col])

run(sys.argv[1],sys.argv[2],int(sys.argv[3]))

