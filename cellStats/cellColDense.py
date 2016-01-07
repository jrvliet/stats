import numpy as np
import matplotlib.pyplot as plt

galIDList = ['D9o2', 'D9q', 'D9m4a']
ionList = ['HI', 'MgII', 'CIV', 'OVI']
labelList = ['dsSN', 'dsALL_1', 'dwALL_8']

fileLoc = '/Users/jacob/research/dwarfs/abscells/'
kpctocm = 3.086e21

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10.5, 10.5))
axes = [ax1, ax2, ax3, ax4]

for galID, lab in zip(galIDList, labelList):

    for ion, ax in zip(ionList,axes):

        column = []

        filename = '{0:s}.{1:s}.bulk_abscells_extend.dat'.format(galID, ion)
        f = open(fileLoc+filename)
        f.readline()
        for line in f:
            l = line.split()
            size = float(l[9])          # In kpc
            ionDense = float(l[14])

            sizecm = size*kpctocm

            column.append( np.log10(sizecm*ionDense) )

        f.close()

        ax.hist(column, bins=20, log=True, histtype='step', label=lab)

for ax, ion in zip(axes,ionList):
    ax.set_xlabel('{0:s} Column Density'.format(ion))
    ax.set_ylabel('Frequency')
    ax.set_xlim([8,22])
    ax.set_ylim([0.1,1e6])
    ax.legend(frameon=False)
    


plt.savefig('colDenseHist.pdf')
