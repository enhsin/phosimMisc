import sys
import numpy as np

"""
http://en.wikipedia.org/wiki/Sellmeier_equation
http://www.researchgate.net/profile/Laurent_Pilon/publication/5819871_Optical_constants_of_silica_glass_from_extreme_ultraviolet_to_far_infrared_at_near_room_temperature/links/0f31752edaeec4acb3000000.pdf

Optical constants of silica glass from extreme ultraviolet
to far infrared at near room temperature
"""
def refractive(l):
    l2 = l*l
    t1 = 0.6961663*l2/(l2-0.0684043**2)
    t2 = 0.4079426*l2/(l2-0.1162414**2)
    t3 = 0.8974794*l2/(l2-9.896161**2)
    n = (1+t1+t2+t3)**0.5
    return n

def makeTable():
    out=open('silica_dispersion.txt','w')
    out.write("# media file for fused silica (corning)\n")
    out.write("# column 1 is wavelength (microns) and column 2 is index of refraction\n")
    for l in np.linspace(0.3, 1.2, 901):
        out.write('    %.5f    %.9f\n' % (l,refractive(l)))
    out.close()

#n=refractive(float(sys.argv[1]))
#print '%.10f' % n
makeTable()
