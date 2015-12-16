
# Plots the phsyical location of cells that give rise to set
# combinations of absorptions

import numpy as np
import matplotlib.pyplot as plt
import tables as tb
import sys


def binup(x, y, numbins):
    H, xedges, yedges = np.histogram2d( x, y, bins=numbins)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    H = np.ma.masked_where(H==0, H)
    H = np.log10(H)
    H = np.rot90(H)
    H = np.flipud(H)
    return H, xedges, yedges


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


numbins = 50
topRow = [0, 4, 8, 12]
midRow = [1, 5, 9, 13]
botRow = [2, 6, 10, 14]
radialIndex = [3, 7, 11, 15]  # Indicies of subplots that are 1D histograms
ion_list = ['HI', 'MgII', 'CIV', 'OVI']
#ion_colls = [HI, MgII, CIV, OVI]
gasfile = 'cellMaster.dat'
h5filename = gasfile.replace('.dat', '.h5')
h5file = tb.open_file(h5filename, mode='r', title='Dwarf Absorbing Cells')
table = h5file.root.cell.readout


for i in range(0,len(ion_list)):
    
    fig,((p11,p12,p13,p14),
         (p21,p22,p23,p24),
         (p31,p32,p33,p34),
         (p41,p42,p43,p44)) = plt.subplots(4,4,figsize=(12,12))
    
    plot_list = [p11, p21, p31, p41,
                 p12, p22, p32, p42,
                 p13, p23, p33, p43,
                 p14, p24, p34, p44]


    ion1 = ion_list[i]
    xed, yed, histos = [], [], []
    
    for j in range(0,len(ion_list)):
        ion2 = ion_list[j]

        # Create the condition
        condition = '({0:s}==1) & ({1:s}==1)'.format(ion1, ion2)
        xloc = [x['xPos'] for x in table.where(condition) ]
        yloc = [x['yPos'] for x in table.where(condition) ]
        zloc = [x['zPos'] for x in table.where(condition) ]
        
        print 'Query finished for {0:s} and {1:s}'.format(ion1, ion2)

        # Bin up x & y
        H, xedges, yedges = binup(xloc, yloc, numbins)
        xed.append(xedges)
        yed.append(yedges)
        histos.append(H)
        
        # Bin up y & z
        H, xedges, yedges = binup(zloc, yloc, numbins)
        xed.append(xedges)
        yed.append(yedges)
        histos.append(H)
        
        # Bin up x & z
        H, xedges, yedges = binup(xloc, zloc, numbins)
        xed.append(xedges)
        yed.append(yedges)
        histos.append(H)
                          
        # Bin up r
        rloc = []
        for k in range(0,len(xloc)):
            x = xloc[k]
            y = yloc[k]
            z = zloc[k]
            r = np.sqrt( x*x + y*y + z*z)
            rloc.append(r)
        H, edges = np.histogram( rloc, bins=numbins)
        xed.append(edges)
        yed.append(edges)
        histos.append(H)
                          
        print 'Binning done'
        
    # Plot the data
    for j in range(0, len(plot_list)):
        print 'Plotting started for ', ion1
        ax = plot_list[j]
        H = histos[j]
        xedges = xed[j]
        yedges = yed[j]

        if j not in radialIndex:
            mesh = ax.pcolormesh(xedges, yedges, H)
        else:
            center = (xedges[:-1] + xedges[1:]) / 2
            width = 0.7 * (xedges[1] - xedges[0])
            ax.bar( center, H, align='center', width=width)

        # Labels:
        if j == topRow[0]:
            ax.set_title( ion1 + ' and ' + ion_list[j/4] )
            ax.set_xlabel( 'x' )
            ax.set_ylabel( 'y' ) 
        if j==topRow[1]:
            ax.set_title( ion1 + ' and ' + ion_list[j/4] )
        if j==topRow[2]:
            ax.set_title( ion1 + ' and ' + ion_list[j/4] )
        if j==topRow[3]:
            ax.set_title( ion1 + ' and ' + ion_list[j/4] )
        if j == midRow[0]:
            ax.set_xlabel( 'z' )
            ax.set_ylabel( 'y' )
        if j == botRow[0]:
            ax.set_xlabel( 'x' )
            ax.set_ylabel( 'y' )
        if j in radialIndex:
            ax.set_xlabel( 'r' )
            ax.set_xlim([0,300])
            if j == radialIndex[0]:
                ax.set_ylabel( 'Number' )
        else:
            ax.set_xlim([-200,200])
            ax.set_ylim([-200,200])

    s = ion1 + 'SpatialLoc.png'
    plt.savefig(s, bbox_inches='tight')
    print '{0:s} saved\n'.format(s)



sys.exit()




f = open('cellMaster.dat')
f.readline()

x, y, z, lognH, logT, massF, cellSize, hi, mgii, civ, ovi = [], [], [], [], [], [], [], [], [], [], []

print 'Read in begin'
linecount = 0
maxlinecount = 1e9
for line in f:
    linecount += 1
#    if linecount < 100:
    l = line.split()
    x.append( float( l[3] ) )
    y.append( float( l[4] ) )
    z.append( float( l[5] ) )
    lognH.append( float( l[6] ) )
    logT.append( float( l[7] ) )
    massF.append( float( l[8] ) )
    cellSize.append( float( l[9] ) )
    hi.append( float( l[10] ) )
    mgii.append( float( l[11] ) )
    civ.append( float( l[12] ) )
    ovi.append( float( l[13] ) )
#    else:
#        break
print len(civ)
f.close()
absorption = [hi, mgii, civ, ovi]
print 'Read in finished'

for i in range(0,len(ion_list)):

    
    ion1 = ion_list[i]
    print ion1
    for j in range(0, len(ion_list)):

        ion2 = ion_list[j]

        xloc = []
        yloc = []
        zloc = []
        rloc = []

        for k in range(0,len(hi)):

            if absorption[i][k]==1 and absorption[j][k]==1:
                xloc.append(x[k])
                yloc.append(y[k])
                zloc.append(z[k])


    
 
