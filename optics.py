class Surface(object):
    def __init__(self,name,curvR,z,rout,rin,conic,a):
        self.name=name
        self.curvR=curvR
        self.z=z
        self.rout=rout
        self.rin=rin
        self.conic=conic
        self.a=a
    def height(self,r):
        if self.curvR==0:
            z=self.z
        else:
            z=self.z-r**2/self.curvR/(1.0+(1.0-(self.conic+1.0)*r**2/self.curvR**2)**0.5)
            for i in range(8):
                z-=self.a[i]*1e3*r**(i+3)
        return z


def readOptics(opticsFile):
    surface=[]
    height=0.0
    for line in open(opticsFile).readlines():
        if "#" in line:
            continue
        a=map(float,line.split()[2:15])
        height+=a[1]
        surfaceName=line.split()[0]
        if surfaceName=='none':
            continue
        radiusofcurv=-a[0]
        outerRadius=a[2]
        innerRadius=a[3]
        conic=a[4]
        surface.append(Surface(surfaceName,radiusofcurv,height,outerRadius,innerRadius,conic,a[5:13]))
    return surface

