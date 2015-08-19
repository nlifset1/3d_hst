print "start"
from astropy.io import ascii
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pylab
import math
import random
import os
from astropy.table import Table
from astropy.cosmology import WMAP9 as cosmo
os.chdir('C:\\3d_hst')

#bring in the data#
data = ascii.read("values_color.dat")

lst_errorsrsr = []
lst_errorsbsr = []
lst_errorsrsb = []
lst_errorsbsb = []

lst_radius = 10**np.linspace(1.2,3.6,13)

lst_galr = []
lst_galb = []

print 'binning'
for gal in data:
    uv = gal['uv']
    vj = gal['vj']
    if (((vj < 0.92) & (uv > 1.3)) | ((vj > 0.8) & (vj < 1.6) & (uv > (0.88*vj +0.49)))):
        lst_galr.append([gal['id'],gal['field']])
    elif (((vj < 0.92) & (uv < 1.3)) | ((vj > 0.8) & (vj < 1.6) & (uv < (0.88*vj +0.49))) | (vj>1.5)):
        lst_galb.append([gal['id'],gal['field']])
    
for i in range(len(lst_radius)):
    print '%s' % (i)
    lst_randrsr = []
    lst_randrsb = []
    for gal in lst_galr:
        gal_info = data[(data['id'] == gal[0]) & (data['field'] == gal[1])]
        lst_randrsr.append(gal_info['randrr{}'.format(i+1)])
        lst_randrsb.append(gal_info['randbr{}'.format(i+1)])
    lst_errorsrsr.append(np.std(lst_randrsr))
    lst_errorsrsb.append(np.std(lst_randrsb))

    lst_randbsr = []
    lst_randbsb = []
    for gal in lst_galb:
        gal_info = data[(data['id'] == gal[0]) & (data['field'] == gal[1])]
        lst_randbsr.append(gal_info['randrr{}'.format(i+1)])
        lst_randbsb.append(gal_info['randbr{}'.format(i+1)])
    lst_errorsbsr.append(np.std(lst_randbsr))
    lst_errorsbsb.append(np.std(lst_randbsb))


data = Table([lst_errorsrsr,lst_errorsrsb,lst_errorsbsr,lst_errorsbsb], names=['stdrsr','stdrsb','std1bsr','stdbsb'])
ascii.write(data, 'color_errors.dat')



