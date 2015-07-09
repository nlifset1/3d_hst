print "start"
from astropy.io import ascii
import numpy as np
from astropy.io import fits
import astropy.constants
import matplotlib.pyplot as plt
import pylab
import math
import random
from scipy import spatial

#bring in the data#
data = ascii.read("aegis_3dhst.v4.1.cat")
data_fast = ascii.read("aegis_3dhst.v4.1.fout")
data_z = ascii.read("aegis_3dhst.v4.0.sfr")

#flag out the bad stuff#
idx, = np.where((data["use_phot"] == 1.0) & (data["star_flag"] == 0.0))
data_fast_flag = data_fast[idx]
data_flag = data[idx]
data_z_flag = data_z[idx]

idx2, = np.where(data_fast_flag["lsfr"] != -99)
data_fast_flagged = data_fast_flag[idx2]
data_flagged = data_flag[idx2]
data_z_flagged = data_z_flag[idx2]

#creating a function for finding number of galaxies within a radius R#
def Counts(gal_id, R):
    #creating a redshift range about the chosen galaxy#
    z_un = data_z_flagged[(data_z_flagged['id'] == gal_id)]
    z_und = z_un['z']
    z = z_und[0]
    #making a list of galaxies in that range#
    lst_gal = []
    for gal in data_fast_flagged:
        if ((gal['z'] >= (z - 0.03)) and (gal['z'] <= (z + 0.03))):
            if gal['id'] != gal_id:
                lst_gal.append(gal['id'])
    #making lists of the euclidean coordinates of every galaxy in the range#
    lst_x = []
    lst_y = []
    lst_z = []
    for gal in lst_gal:
        #pulling the necessary info#
        position_info = data_flagged[(data_flagged['id'] == gal)]
        z_info = data_z_flagged[(data_z_flagged['id'] == gal)]

        ra_ = (position_info['ra'])*0.0174532925
        ra = ra_[0]
        dec_ = (position_info['dec'])*0.0174532925
        dec = dec_[0]
        dist_ = (z_info['z'])*4285
        dist = dist_[0]
        #calculating the euclidean coordinates#
        x = (math.cos(ra))*(math.cos(dec))*(dist)
        y = (math.sin(ra))*(math.cos(dec))*(dist)
        z = (math.sin(ra))*(dist)

        #putting those coordinates into the lists#
        lst_x.append(x)
        lst_y.append(y)
        lst_z.append(z)
    #making a KD Tree of the euclidean info#
    points = zip(lst_x, lst_y, lst_z)
    tree = spatial.KDTree(points)
    #getting euclidean coordinates of the specific galaxy#
    p1 = data_flagged[(data_flagged['id'] == gal_id)]
    ra1_ = (p1['ra'])*0.0174532925
    ra1 = ra1_[0]
    dec1_ = (p1['dec'])*0.0174532925
    dec1 = dec1_[0]
    dist1_ = (z_un['z'])*4282.7494
    dist1 = dist1_[0]
    x1 = (math.cos(ra1))*(math.cos(dec1))*(dist1)
    y1 = (math.sin(ra1))*(math.cos(dec1))*(dist1)
    z1 = (math.sin(ra1))*(dist1)
    #finding number of galaxies within radius R of the specific galaxy#
    within_lst = tree.query_ball_point([x1,y1,z1], R)
    within = len(within_lst)
    return within

#creating a function similar to Counts but for random base line#
def rand_counts(R):
    #picking random location for galaxy number density#
    ra1 = random.uniform(3.746000, 3.756821)
    dec1 = random.uniform(0.920312, 0.925897)
    z = random.uniform(0.4, 2.0)
    #making a list of galaxies within a specific redshift range of random point#
    lst_gal = []
    for gal in data_fast_flagged:
        if ((gal['z'] >= (z - 0.03)) and (gal['z'] <= (z + 0.03))):
                lst_gal.append(gal['id'])
    #making lists of the euclidean coordinates of galaxies in that range#
    x_lst = []
    y_lst = []
    z_lst = []
    for gal in lst_gal:
        #pulling the necessary info#
        position_info = data_flagged[(data_flagged['id'] == gal)]
        z_info = data_z_flagged[(data_z_flagged['id'] == gal)]

        ra_ = (position_info['ra'])*0.0174532925
        ra = ra_[0]
        dec_ = (position_info['dec'])*0.0174532925
        dec = dec_[0]
        dist_ = (z_info['z'])*4285
        dist = dist_[0]
        #calculating the euclidean coordinates#
        x = (math.cos(ra))*(math.cos(dec))*(dist)
        y = (math.sin(ra))*(math.cos(dec))*(dist)
        z = (math.sin(ra))*(dist)

        #putting those coordinates into the lists#
        x_lst.append(x)
        y_lst.append(y)
        z_lst.append(z)
    #making a KD Tree of the info#
    points = zip(x_lst, y_lst, z_lst)
    tree = spatial.KDTree(points)
    #getting euclidean coordinates of specific galaxy#
    dist1 = (z)*4282.7494
    x1 = (math.cos(ra1))*(math.cos(dec1))*(dist1)
    y1 = (math.sin(ra1))*(math.cos(dec1))*(dist1)
    z1 = (math.sin(ra1))*(dist1)
    #calculating number of galaxies within radius R#
    within_lst = tree.query_ball_point([x1,y1,z1], R)
    within = len(within_lst)
    return within

#getting list of galaxies that avoid edges by 0.05 degrees#
lst_gal_edged =[]
for gal in data_flagged:
    if ((gal['ra'] < 215.25) and (gal['ra'] > 214.63)):
        if ((gal['dec'] < 53.05) and (gal['dec'] > 52.73)):
            lst_gal_edged.append(gal['id'])


#getting a list of galaxies with lmass >= 11.0 and within the redshift range of 0.4<z<2#
lst_gal_massed = []
for gal in lst_gal_edged:
    gal_info = data_fast_flagged[(data_z_flagged['id'] == gal)]
    if ((gal_info['z'] >= 0.4) & (gal_info['z'] <= 2.0)):
        if (gal_info['lmass'] >= 11.0):
            lst_gal_massed.append(gal)

#getting a list of all galaxies that avoid edges and aren't massive#
lst_gal_small = lst_gal_edged
for gal in lst_gal_massed:
    lst_gal_small.remove(gal)

#making lists for the plots of mass vs density, R is in MPC#
R = 20
rand_counted = 0
for i in range(len(lst_gal_massed)):
    rand_counted += rand_counts(R)
rand_counted = float(rand_counted/len(lst_gal_massed))
lst_z =[]
lst_counts =[]
for gal in lst_gal_massed:
    gal_counted = Counts(gal, R)
    lst_counts.append((float(gal_counted) - rand_counted)/((R**2)*math.pi))
    for item in data_z_flagged:
        if item['id'] == gal:
            lst_z.append(item['z'])


#plotting redshift vs density#

pylab.scatter(lst_z, lst_counts)

pylab.suptitle('Redshift vs Log Galaxy Number Density of All Redshifts', fontsize=20)
pylab.xlabel('Redshift', fontsize=16)
pylab.ylabel('Log Galaxy Number Density (mpc^-2)', fontsize=15)
pylab.xlim([0.4,2.0])
pylab.ylim([0.004,0.06])
pylab.yscale('log')


pylab.ion()
pylab.show()

