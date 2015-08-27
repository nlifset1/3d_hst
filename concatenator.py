from astropy.io import ascii
import numpy as np
import os
from astropy.table import Table, vstack
os.chdir('C:\\3d_hst')


#bring in the data#
print 'aegis'
data_aegis = ascii.read('C:\\3d_hst/aegis_3dhst.v4.1.master.RF')
print 'cosmos'
data_cosmos = ascii.read('C:\\3d_hst/cosmos_3dhst.v4.1.master.RF')
print 'goodss'
data_goodss = ascii.read('C:\\3d_hst/goodss_3dhst.v4.1.master.RF')
print 'goodsn'
data_goodsn = ascii.read('C:\\3d_hst/goodsn_3dhst.v4.1.master.RF')
print 'uds'
data_uds = ascii.read('C:\\3d_hst/uds_3dhst.v4.1.master.RF')

#concatenating all the fields into one and writing it to a file#
print 'concatenating'
data_1 = vstack([data_aegis, data_cosmos])
data_2 = vstack([data_goodss, data_goodsn])
data_3 = vstack([data_1, data_2])
data_big = vstack([data_3, data_uds])
ascii.write(data_big, '3dhst_big.dat')
