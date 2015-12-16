

import numpy as np
import matplotlib.pyplot as plt

lab_list = ['dwSN', 'dwALL_1', 'dwALL_8']
galID_list = ['D9o2', 'D9q', 'D9m4a']
ion_list = ['HI', 'MgII', 'CIV', 'OVI']

absLoc = '/home/jacob/research/dwarfs/abscells/'

fig, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2,2,figsize=(10.2,10.2))
axes = [ax11, ax12, ax21, ax22]

for ion in ion_list:
    
    ax = axes[ion_list.index(ion)]
    print ion
    for galID,lab in zip(galID_list, lab_list):

        filename = '{0:s}.{1:s}.bulk_abscells.dat'.format(galID,ion)
        f = open(absLoc+filename)
        f.readline()
        logN = []
        for line in f:
            logN.append(float(line.split()[4]))
        f.close()

        ax.hist(logN, bins=20, log=True, histtype='step', label=lab)

    ax.legend(frameon=False)
    ax.set_xlabel('log N$_{ion}$')
    ax.set_ylabel('Frequency')
    ax.set_title(ion)
fig.tight_layout()
s = 'logNHist.png'
fig.savefig(s)








