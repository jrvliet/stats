#
# Determines the distribution of absorbing cells along a LOS
# Read in the list of absorbing cellIDs from the abscells files
# Get the 

import sys
import numpy as np
import random as rn
import matplotlib.pyplot as plt
import math

def numNeighbors(points, length, numbins):
    
    # Determines the mean number of neighbors within certain distance 
    # of each point. Returns an array of number of neighbors

    numsteps = numbins
    step = length/numsteps
    
    numNeigh = np.zeros(numsteps)

    # Fill the distance limits 
    minDist = np.zeros(numsteps)
    maxDist = np.zeros(numsteps)
    for i in range(numsteps):
        minDist[i] = i*step
        maxDist[i] = (i+1)*step
    
    # Loop over all points in the array
    for i in range(0,len(points)):
        
        # Get the location of the prime point
        primeLoc = points[i]

        # Loop over all other points
        for j in range(i,len(points)):
            # Don't include the original point
            if j!=i:
                # Get the distance between this point and the prime location
                dist = points[j]-primeLoc
                
                # Get the distance bin that this pair belongs in
                for k in range(0,len(minDist)):
                    if dist>minDist[k] and dist<maxDist[k]:
                        numNeigh[k] += 1

    return numNeigh

numbins = 10

ion_list = ['HI', 'MgII', 'CIV', 'OVI']
galID_list = ['D9o2', 'D9q', 'D9m4a']
expn_list1 = '0.900 0.926 0.950 0.975 0.990 1.002'.split()
expn_list2 = '0.901 0.925 0.950 0.976 0.991 1.001'.split()
expn_list3 = '0.900 0.925 0.950 0.975 0.990 1.000'.split()
expn_list = [expn_list1, expn_list2, expn_list3]

# For testing
#ion_list = ['MgII']
#galID_list = ['D9o2']
#expn_list = [expn_list1]

#fig, ((ax11,ax12),(ax21,ax22)) = plt.subplots(2,2,figsize=(10.2,10.2))
#axes = [ax11, ax12, ax21, ax22]

abscellLoc = '/home/jacob/research/dwarfs/abscells/individual/'
gasboxLoc = '/home/jacob/research/dwarfs/gasfiles/'
linesLoc = '/home/jacob/research/dwarfs/lines/'

for ion in ion_list:

    print ion
    ratio = np.zeros(numbins)
    rcount = np.zeros(numbins)
    avelengths = []
    for galID, expn in zip(galID_list, expn_list):
        print '\t', galID
        for a in expn:
            print '\t\t',a
            # Read in lines files
            linesName = '{0:s}/{1:s}_{2:s}_lines.dat'.format(linesLoc,galID,a)
            lines = np.loadtxt(linesName, skiprows=2)

            # Read in gas box
            gasboxName = '{0:s}/{1:s}_GZa{2:s}.txt'.format(gasboxLoc,galID,a)
            gasbox = np.loadtxt(gasboxName, skiprows=2) 
            # x=1, y=2, z=3

            # Read in abscell data
            abscellName = '{0:s}/{1:s}.{2:s}.{3:s}.abs_cells.dat'.format(
                            abscellLoc, galID, a, ion)
            fabs = open(abscellName)
            fabs.readline()
    
            # along is a 2d array
            # each row is a LOS
            # along each row will be the s for each absorbing cell
            along = np.zeros((999,400))
            numalong = np.zeros(999)
            dumcount = 0    
            
            # Loop over all absorbing cells
            for line in fabs:
                dumcount+=1
                losnum = int(line.split()[0])
                cellnum = int(line.split()[2])

                xen = lines[losnum-1][0]
                yen = lines[losnum-1][1]
                zen = lines[losnum-1][2]

                x = gasbox[cellnum-1][1]*1000
                y = gasbox[cellnum-1][2]*1000
                z = gasbox[cellnum-1][3]*1000
            
                # Calculated the distance the cell is along the LOS
                s = np.sqrt( pow(x-xen,2)+pow(y-yen,2)+pow(z-zen,2) )
#                print losnum, cellnum, xen, yen, zen, x, y, z, s
                
                # Add to the along array
                currentInd = numalong[losnum-1]
                along[losnum-1][currentInd] = s
                numalong[losnum-1] += 1

#                if dumcount>3:
#                    sys.exit()
            
            # Get the length of each LOS
#            length = np.zeros(999)
            length = []
            for i in range(0,len(lines[:,0])):
                xen = lines[i][0]
                yen = lines[i][1]
                zen = lines[i][2]
                
                xex = lines[i][3]
                yex = lines[i][4]
                zex = lines[i][5]

                l = np.sqrt(pow(xen-xex,2) + pow(yen-yex,2) + pow(zen-zex,2) ) 
#                length[i] = l
                length.append(l)
#            print 'Length mean: {0:f} max: {1:f} min: {2:f}'.format(np.mean(length),
#                    np.max(length), np.min(length))
            
#            print 'Numalong mean: {0:f} max: {1:f} min: {2:f}'.format(np.mean(numalong),
#                    np.max(numalong), np.min(numalong))
            
            # Get the average distance between pairs of 
            # absorbing cells for a LOS
            # Loop over LOS
            ratiof = open('ratio.out', 'w')
            for i in range(0,len(numalong)):
            
                numberofCells = int(numalong[i])
                if numberofCells>0:
                    l = length[i]
                    cellLocs = along[i]
                    dcell = []
                    drand = []
                    randLoc = []

                    pointsCorrelation = numNeighbors(cellLocs, l, numbins)


#                    print cellLocs
#                    sys.exit()
                    # Get the distance between pairs of cells
#                    for j in range(0,numberofCells):
#                        for k in range(j+1, numberofCells):
#                            distance = abs( cellLocs[j]-cellLocs[k] )
#                            dcell.append(distance)
#                            print distance

                    # Get the distance between pairs of random cells
                    for j in range(0,numberofCells):
                        rL = rn.random()*l
                        randLoc.append(rL)
                    randCorrelation = numNeighbors(randLoc, l, numbins)

#                    for j in range(0,numberofCells):
#                        for k in range(j+1, numberofCells):
#                            distance = abs( randLoc[j]-randLoc[k] )
#                            drand.append(distance)
                    
#                    meandcell = np.mean(dcell)
#                    meandrand = np.mean(drand)
#                    print meandcell, meandrand

                    avelengths.append(l)
                    for j in range(0,len(ratio)):
                        try:
                            r = pointsCorrelation[j]/randCorrelation[j]
                            if not math.isnan(r) and not math.isinf(r):
                                ratiof.write('{0:f}\t{1:f}\t{2:f}\n'.format(pointsCorrelation[j], randCorrelation[j],r))
                                ratio[j] += r
                            
                        except RuntimeWarning:
                            r = 0
                        
                        rcount[j] += 1
                

    for i in range(0,len(ratio)):
        ratio[i] = ratio[i]/rcount[i]

    print ratio
    averageL = sum(avelengths)/len(avelengths)
    print 'Average LOS lenght: {0:f}'.format(averageL)
    binlims = []
    binsteps = averageL/numbins
    for i in range(0,numbins):
        binlims.append(i*binsteps)

    print np.mean(ratio)
    print np.max(ratio)
    print np.mean(ratio) 
    # Plot histogram of ratios
#    ax = axes[ion_list.index(ion)]
#    ax.hist(ratio, bins=20, log=True, histtype='step')
#    ax.plot(binlims, ratio, 'kx')
#    ax.set_title(ion)
#    ax.set_xlim([0,600])
#    ax.set_xlabel('Distance [kpc]')
#    ax.set_ylabel('Ratio')

    plt.plot(binlims, ratio, label=ion)
    plt.xlabel('Distance [kpc]')
    plt.ylabel('Ratio')
    

#plt.tight_layout()
s = 'clumping.pdf'
plt.savefig(s)




            

    
