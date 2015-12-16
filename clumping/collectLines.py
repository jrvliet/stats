#
# Collect lines.dat files
#

import numpy as np
import subprocess as sp

ion_list = ['HI', 'MgII', 'CIV', 'OVI']
galID_list = ['D9o2', 'D9q', 'D9m4a']
expn_list1 = '0.900 0.926 0.950 0.975 0.990 1.002'.split()
expn_list2 = '0.901 0.925 0.950 0.976 0.991 1.001'.split()
expn_list3 = '0.900 0.925 0.950 0.975 0.990 1.000'.split()
expn_list = [expn_list1, expn_list2, expn_list3]


baseLoc = '/home/matrix3/jrvander/sebass_gals/dwarfs/'
newLoc = baseLoc+'lines/'

for galID, expn in zip(galID_list, expn_list):

    for a in expn:
    
        loc = '{0:s}/{1:s}_outputs/a{2:s}/'.format(baseLoc,galID,a)
        oldname = 'lines.dat'
        newname = '{0:s}_{1:s}_lines.dat'.format(galID, a)
        command = 'cp {0:s}{1:s} {2:s}{4}'.format(loc,oldname,newLoc,newname)
        
        sp.call(command, shell=True)




   
