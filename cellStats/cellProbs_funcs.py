#
# Functions to help cellProbs.py 
# Mostly with writing to file
#

import numpy as np

def printCellProperties( fo, ion, density, temp, metal, size):

    # Write a block to outFile with a summaries of the 
    # properties passed in

    denseMean = np.log10(np.mean(density))
    tempMean  = np.log10(np.mean(temp))
    metalMean = np.log10(np.mean(metal))
    sizeMean  = np.log10(np.mean(size))

    denseMin = np.log10(np.min(density))
    tempMin  = np.log10(np.min(temp))
    metalMin = np.log10(np.min(metal))
    sizeMin  = np.log10(np.min(size))

    denseMax = np.log10(np.max(density))
    tempMax  = np.log10(np.max(temp))
    metalMax = np.log10(np.max(metal))
    sizeMax  = np.log10(np.max(size))

    denseMed = np.log10(np.median(density))
    tempMed  = np.log10(np.median(temp))
    metalMed = np.log10(np.median(metal))
    sizeMed  = np.log10(np.median(size))

    denseStd = np.log10(np.std(density))
    tempStd  = np.log10(np.std(temp))
    metalStd = np.log10(np.std(metal))
    sizeStd  = np.log10(np.std(size))

    s = 'Typical {0:s} cell properties\n'.format(ion)
    fo.write(s)
    s = 'Property \t Min \t Max \t Mean \t Med \t Std Dev \n'
    fo.write(s)
    s = 'Density \t {0:.2f} \t {1:.2f} \t {2:.2f} \t {3:.2f} \t {4:.2f} \n'.format(denseMin, denseMax, denseMean, denseMed, denseStd)
    fo.write(s)
    s = 'Temperature \t {0:.2f} \t {1:.2f} \t {2:.2f} \t {3:.2f} \t {4:.2f} \n'.format(tempMin, tempMax, tempMean, tempMed, tempStd)

    s = 'Metalicity \t {0:.2f} \t {1:.2f} \t {2:.2f} \t {3:.2f} \t {4:.2f} \n'.format(metalMin, metalMax, metalMean, metalMed, metalStd)
    fo.write(s)
    s = 'Cell Size \t {0:.2f} \t {1:.2f} \t {2:.2f} \t {3:.2f} \t {4:.2f} \n'.format(sizeMin, sizeMax, sizeMean, sizeMed, sizeStd)
    fo.write(s)
    s = 'Fraction of MgII absorbing cells that do not have HI absorption: {0:.2%}\n'.format(mgFrac)
    fo.write(s)




def findCellProps( datfile, ion):

    # Find the properties of all cells that do not
    # have significant absorpiton in the ion passed in

    ionList = ['HI', 'MgII', 'CIV', 'OVI']
    target = ionList.index(ion) + 1

    f = open('dat.files')  # A file that is a list of file names
    ion1_nH, ion1_t, ion1_z, ion1_l = [], [], [], []
    ion2_nH, ion2_t, ion2_z, ion2_l = [], [], [], []
    ion3_nH, ion3_t, ion3_z, ion3_l = [], [], [], []
    ion4_nH, ion4_t, ion4_z, ion4_l = [], [], [], []
    ion1Count = 0.
    ion2Count = 0.
    ion3Count = 0.
    count = [ion1Count, ion2Count, ion3Count]
    nH = [ion1_nH, ion2_nH, ion3_nH, ion4_nH]
    t = [ion1_t, ion2_t, ion3_t, ion4_t]
    z = [ion1_z, ion2_z, ion3_z, ion4_z]
    l = [ion1_l, ion2_l, ion3_l, ion4_l]

    for fN in f:
        ftmp = open(fN.strip())
        ftmp.readline()

        for line in ftmp:
            l = line.split()
            hiFlag = int(l[8])
            mgiiFlag = int(l[9])
            civFlag = int(l[10])
            oviFlag = int(l[11])
            
            flags = [l[8], l[9], l[10], l[11]]
        
            for i in range(0, 4):
                if i!=target:
                    if flags[target]==0 and flags[i]==1:
                        nH[i].append( pow(10, float(l[4])))     
                        t[i].append( pow(10, float(l[5])))     
                        z[i].append( pow(10, float(l[6])))     
                        l[i].append( pow(10, float(l[7])))     
                        count[i] += 1
        ftmp.close()
    f.close()

    ind = [0, 1, 2, 3]
    ind.pop(target)
    
    ion1 = [nH[ind[0]], t[ind[0]], z[ind[0]], l[ind[0]]]
    ion2 = [nH[ind[1]], t[ind[1]], z[ind[1]], l[ind[1]]]
    ion3 = [nH[ind[2]], t[ind[2]], z[ind[2]], l[ind[2]]]

    return ion1, ion2, ion3

