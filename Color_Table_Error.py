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

#initializing lists#
#a little help: lst_errors [color of central] s [color of satellites] [redshift bin]#
lst_errorsrsr1 = []
lst_errorsbsr1 = []
lst_errorsrsb1 = []
lst_errorsbsb1 = []

lst_errorsrsr2 = []
lst_errorsbsr2 = []
lst_errorsrsb2 = []
lst_errorsbsb2 = []

lst_errorsrsr3 = []
lst_errorsbsr3 = []
lst_errorsrsb3 = []
lst_errorsbsb3 = []

lst_errorsrsr4 = []
lst_errorsbsr4 = []
lst_errorsrsb4 = []
lst_errorsbsb4 = []

lst_radius = 10**np.linspace(1.2,3.6,13)

lst_galr = []
lst_galb = []

#separating the galaxies into quiescent(r) and star-forming(b) lists#
print 'binning'
for gal in data:
    uv = gal['uv']
    vj = gal['vj']
    if (((vj < 0.92) & (uv > 1.3)) | ((vj > 0.8) & (vj < 1.6) & (uv > (0.88*vj +0.49)))):
        lst_galr.append([gal['id'],gal['field']])
    elif (((vj < 0.92) & (uv < 1.3)) | ((vj > 0.8) & (vj < 1.6) & (uv < (0.88*vj +0.49))) | (vj>1.5)):
        lst_galb.append([gal['id'],gal['field']])

#splitting the central galaxies into more bins based on redshift#
#low number is low redshift#
lst_galr1 = []
lst_galb1 = []
lst_galr2 = []
lst_galb2 = []
lst_galr3 = []
lst_galb3 = []
lst_galr4 = []
lst_galb4 = []
for gal in lst_galr:
    gal_info = data[(data['id'] == gal[0]) & (data['field'] == gal[1])]
    if gal_info['z'] < 1.0:
        lst_galr1.append(gal)
    elif gal_info['z'] < 1.5:
        lst_galr2.append(gal)
    elif gal_info['z'] < 2.0:
        lst_galr3.append(gal)
    elif gal_info['z'] < 2.5:
        lst_galr4.append(gal)

for gal in lst_galb:
    gal_info = data[(data['id'] == gal[0]) & (data['field'] == gal[1])]
    if gal_info['z'] < 1.0:
        lst_galb1.append(gal)
    elif gal_info['z'] < 1.5:
        lst_galb2.append(gal)
    elif gal_info['z'] < 2.0:
        lst_galb3.append(gal)
    elif gal_info['z'] < 2.5:
        lst_galb4.append(gal)

#for each radius and then each galaxy, pulling data on the random densities, and calculating std for error bars#
for i in range(len(lst_radius)):
    print '%s' % (i)
    lst_randrsr1 = []
    lst_randrsb1 = []
    for gal in lst_galr1:
        gal_info = data[(data['id'] == gal[0]) & (data['field'] == gal[1])]
        lst_randrsr1.append(gal_info['randrr{}'.format(i+1)])
        lst_randrsb1.append(gal_info['randbr{}'.format(i+1)])
    lst_errorsrsr1.append(np.std(lst_randrsr1))
    lst_errorsrsb1.append(np.std(lst_randrsb1))

    lst_randbsr1 = []
    lst_randbsb1 = []
    for gal in lst_galb1:
        gal_info = data[(data['id'] == gal[0]) & (data['field'] == gal[1])]
        lst_randbsr1.append(gal_info['randrr{}'.format(i+1)])
        lst_randbsb1.append(gal_info['randbr{}'.format(i+1)])
    lst_errorsbsr1.append(np.std(lst_randbsr1))
    lst_errorsbsb1.append(np.std(lst_randbsb1))

    lst_randrsr2 = []
    lst_randrsb2 = []
    for gal in lst_galr2:
        gal_info = data[(data['id'] == gal[0]) & (data['field'] == gal[1])]
        lst_randrsr2.append(gal_info['randrr{}'.format(i+1)])
        lst_randrsb2.append(gal_info['randbr{}'.format(i+1)])
    lst_errorsrsr2.append(np.std(lst_randrsr2))
    lst_errorsrsb2.append(np.std(lst_randrsb2))

    lst_randbsr2 = []
    lst_randbsb2 = []
    for gal in lst_galb2:
        gal_info = data[(data['id'] == gal[0]) & (data['field'] == gal[1])]
        lst_randbsr2.append(gal_info['randrr{}'.format(i+1)])
        lst_randbsb2.append(gal_info['randbr{}'.format(i+1)])
    lst_errorsbsr2.append(np.std(lst_randbsr2))
    lst_errorsbsb2.append(np.std(lst_randbsb2))

    lst_randrsr3 = []
    lst_randrsb3 = []
    for gal in lst_galr3:
        gal_info = data[(data['id'] == gal[0]) & (data['field'] == gal[1])]
        lst_randrsr3.append(gal_info['randrr{}'.format(i+1)])
        lst_randrsb3.append(gal_info['randbr{}'.format(i+1)])
    lst_errorsrsr3.append(np.std(lst_randrsr3))
    lst_errorsrsb3.append(np.std(lst_randrsb3))

    lst_randbsr3 = []
    lst_randbsb3 = []
    for gal in lst_galb3:
        gal_info = data[(data['id'] == gal[0]) & (data['field'] == gal[1])]
        lst_randbsr3.append(gal_info['randrr{}'.format(i+1)])
        lst_randbsb3.append(gal_info['randbr{}'.format(i+1)])
    lst_errorsbsr3.append(np.std(lst_randbsr3))
    lst_errorsbsb3.append(np.std(lst_randbsb3))

    lst_randrsr4 = []
    lst_randrsb4 = []
    for gal in lst_galr4:
        gal_info = data[(data['id'] == gal[0]) & (data['field'] == gal[1])]
        lst_randrsr4.append(gal_info['randrr{}'.format(i+1)])
        lst_randrsb4.append(gal_info['randbr{}'.format(i+1)])
    lst_errorsrsr4.append(np.std(lst_randrsr4))
    lst_errorsrsb4.append(np.std(lst_randrsb4))

    lst_randbsr4 = []
    lst_randbsb4 = []
    for gal in lst_galb4:
        gal_info = data[(data['id'] == gal[0]) & (data['field'] == gal[1])]
        lst_randbsr4.append(gal_info['randrr{}'.format(i+1)])
        lst_randbsb4.append(gal_info['randbr{}'.format(i+1)])
    lst_errorsbsr4.append(np.std(lst_randbsr4))
    lst_errorsbsb4.append(np.std(lst_randbsb4))

#writing the table#
data = Table([lst_errorsrsr1,lst_errorsrsb1,lst_errorsbsr1,lst_errorsbsb1,lst_errorsrsr2,lst_errorsrsb2,lst_errorsbsr2,lst_errorsbsb2,lst_errorsrsr3,lst_errorsrsb3,lst_errorsbsr3,lst_errorsbsb3,lst_errorsrsr4,lst_errorsrsb4,lst_errorsbsr4,lst_errorsbsb4], names=['stdrsr1','stdrsb1','std1bsr1','stdbsb1','stdrsr2','stdrsb2','std1bsr2','stdbsb2','stdrsr3','stdrsb3','std1bsr3','stdbsb3','stdrsr4','stdrsb4','std1bsr4','stdbsb4'])
ascii.write(data, 'color_errors.dat')



