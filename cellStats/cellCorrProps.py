
""" Deterimines the typical properites of cells that
    have some absorption properites and not others.
    For example, the typical cell size of cells that 
    have CIV but not MgII absorption
    """


import numpy as np
import matplotlib.pyplot as plt
import tables as tb


def plotColumn(data, column, ion1, ion2, mode):

    """ Plots a column of plots
        From top to bottom: r, L, n, T, Z
        """
    
    numbins = 50
    labels = ['r [kpc]', 'Cell size [kpc]', 'log nH', 'log T [K]', 'SNII mass fraction']
    xmins = [0, 0, -9, 2, -14]
    xmaxs = [300, 30, 2, 8, -2]     


    for i in range(0,len(column)):
        ax = column[i]
        x = data[i]
        ax.hist(x, bins=numbins, range=(xmins[i], xmaxs[i]), histtype='step', log=True)
        ax.set_xlabel(labels[i])
        ax.set_xlim( [ xmins[i], xmaxs[i] ] )
#        ax.set_ylabel('Count')
    
    ax = column[0]
    if mode=='both':
        # Both ions are present
        ax.set_title('{0:s} and {1:s}'.format(ion1, ion2))
    elif mode=='1':
        # Only ion1 is present
        ax.set_title('{0:s} without {1:s}'.format(ion1, ion2))
    elif mode=='2':
        # Only ion2 is present
        ax.set_title('{0:s} without {1:s}'.format(ion2, ion1))
        
        



ions = ['HI', 'MgII', 'CIV', 'OVI']

hi, mgii, civ, ovi = [], [], [], []
x, y, z, n, t, metal, size = [], [], [], [], [], [], []

# Define the columns for the gas box
class cells(tb.IsDescription):
    galID     = tb.StringCol(8) # 8-character string
    expn      = tb.StringCol(8)
    cellID    = tb.UInt32Col()   # Unsigned integer
    xPos      = tb.Float32Col()
    yPos      = tb.Float32Col()
    zPos      = tb.Float32Col()
    lognH     = tb.Float32Col()
    logT      = tb.Float32Col()
    SNII      = tb.Float32Col()
    cellsize  = tb.Float32Col() # float (single-precision) (64 bits)
    HI        = tb.UInt32Col()
    MgII      = tb.UInt32Col()
    CIV       = tb.UInt32Col()
    OVI       = tb.UInt32Col()



filename = 'cellMaster.h5'
h5file = tb.open_file(filename, mode='r', title='Dwarf Absorbing Cells')
table = h5file.root.cell.readout

for i in range(0,len(ions)):
    ion1 = ions[i]
    for j in range(i+1,len(ions)):

        ion2 = ions[j]

        print '{0:s} and {1:s}'.format(ion1, ion2)

        fig, ((p11, p12, p13),
              (p21, p22, p23),
              (p31, p32, p33),
              (p41, p42, p43),
              (p51, p52, p53)) = plt.subplots(5, 3, figsize=(10.5, 10))

        column1 = [p11, p21, p31, p41, p51]
        column2 = [p12, p22, p32, p42, p52]
        column3 = [p13, p23, p33, p43, p53]     

        # Create the condition for both ions being present
        condition = '({0:s}==1) & ({1:s}==1)'.format(ion1, ion2)
        xloc = [x['xPos'] for x in table.where(condition)]
        yloc = [x['yPos'] for x in table.where(condition)]
        zloc = [x['zPos'] for x in table.where(condition)]
        n = [x['lognH'] for x in table.where(condition)]
        t = [x['logT'] for x in table.where(condition)]
        metal = [x['SNII'] for x in table.where(condition)]
        size = [x['cellsize'] for x in table.where(condition)]

        metal = [np.log10(index) for index in metal]

        print 'For condition {0:s}, found {1:d}'.format(condition, len(xloc))
        # Build up galactocentric distance
        rloc = []
        for k in range(0,len(xloc)):
            x = xloc[k]
            y = yloc[k]
            z = zloc[k]
            r = np.sqrt( x*x + y*y + z*z )
            rloc.append(r)

        data = [rloc, size, n, t, metal]
        plotColumn(data, column1, ion1, ion2, 'both')
        print '\tColumn 1 Plotted'

        
        # Search for condition with ion1, but not ion2
        condition = '({0:s}==1) & ({1:s}==0)'.format(ion1, ion2)
        xloc = [x['xPos'] for x in table.where(condition)]
        yloc = [x['yPos'] for x in table.where(condition)]
        zloc = [x['zPos'] for x in table.where(condition)]
        n = [x['lognH'] for x in table.where(condition)]
        t = [x['logT'] for x in table.where(condition)]
        metal = [x['SNII'] for x in table.where(condition)]
        size = [x['cellsize'] for x in table.where(condition)]
        metal = [np.log10(index) for index in metal]

        numFound = len(xloc)
        print 'For condition {0:s}, found {1:d}'.format(condition, len(xloc))
        # Build up galactocentric distance
        rloc = []
        for k in range(0,len(xloc)):
            x = xloc[k]
            y = yloc[k]
            z = zloc[k]
            r = np.sqrt( x*x + y*y + z*z )
            rloc.append(r)

        data = [rloc, size, n, t, metal]
        plotColumn(data, column2, ion1, ion2, '1')
        print '\tColumn 2 Plotted'


        # Search for condition with ion2, but not ion1
        condition = '({0:s}==0) & ({1:s}==1)'.format(ion1, ion2)
        xloc = [x['xPos'] for x in table.where(condition)]
        yloc = [x['yPos'] for x in table.where(condition)]
        zloc = [x['zPos'] for x in table.where(condition)]
        n = [x['lognH'] for x in table.where(condition)]
        t = [x['logT'] for x in table.where(condition)]
        metal = [x['SNII'] for x in table.where(condition)]
        size = [x['cellsize'] for x in table.where(condition)]
        print 'For condition {0:s}, found {1:d}'.format(condition, len(xloc))
        metal = [np.log10(index) for index in metal]

        # Build up galactocentric distance
        rloc = []
        for k in range(0,len(xloc)):
            x = xloc[k]
            y = yloc[k]
            z = zloc[k]
            r = np.sqrt( x*x + y*y + z*z )
            rloc.append(r)

        data = [rloc, size, n, t, metal]
        plotColumn(data, column3, ion1, ion2, '2')
        print '\tColumn 3 Plotted'

                
        
        # Save the plot
        fig.tight_layout()
        s = ion1+'_'+ion2+'_corrProps_log.pdf'
        plt.savefig(s, bbox_inches='tight')
        print '\t{0:s} saved\n'.format(s)



