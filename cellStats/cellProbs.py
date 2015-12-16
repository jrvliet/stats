
#
# Calculates the probabilty of finding an ion absorption
# given the existance of another ion absorption


import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
from cellProbs_funcs import *


# Calculates the probability of finding b given a
def probAgivenB(a, b):

    totalCells = float(len(a))
    probA = float( a.count(1) ) / totalCells
    
    countBoth = 0.
    for i in range(0,len(a)):
        
        if a[i]==1 and b[i]==1:
            countBoth += 1.

    probAandB = countBoth / totalCells
    probAgB = probAandB / probA
    
    return probAgB

hi   = []
mgii = []
civ  = []
ovi  = []

f = open('dat.files')

for fileName in f:

    print fileName.strip()
    fN = fileName.strip()
    
    ftmp = open(fN)
    ftmp.readline()

    for line in ftmp:

        l2 = line.split()
        hiTmp   = float(l2[8])
        mgiiTmp = float(l2[9])
        civTmp  = float(l2[10])
        oviTmp  = float(l2[11])
        
        hi.append(hiTmp)
        mgii.append(mgiiTmp)
        civ.append(civTmp)
        ovi.append(oviTmp)


    ftmp.close()

f.close()


# Read in the list of all probbed cells
f = open('probed.files')
probedCells = []
for fileName in f:

    fN = fileName.strip()
    ftmp = open(fN)
    ftmp.readline()

    for line in ftmp:
    
        cellnum = int(line)
        probedCells.append(cellnum)

    ftmp.close()
f.close()


# Calculate the probabiltiy of finding MgII given that HI is found
totalCells = float(len(probedCells))
#totalCells = float(len(hi))
print ''
print 'Total Number of Cells:          ', totalCells
print 'Number of HI Absorbing Cells:   ', hi.count(1)
print 'Number of MgII Absorbing Cells: ', mgii.count(1)
print 'Number of CIV Absorbing Cells:  ', civ.count(1)
print 'Number of OVI Absorbing Cells:  ', ovi.count(1)
print ''

# Find probability of each ion
probHI   = float(hi.count(1))/totalCells
probMgII = float(mgii.count(1))/totalCells
probCIV  = float(civ.count(1))/totalCells
probOVI  = float(ovi.count(1))/totalCells

# Find the number of occurances of HI and MgII
count = 0.
for i in range(0,len(hi)):

    if hi[i]==1 and mgii[i]==1:
        count+=1

probAandB = count/totalCells

probMgIIgivenHI = probAandB / probHI
print probMgIIgivenHI


probHI = float(hi.count(1))/totalCells 
probMgII = float(mgii.count(1))/totalCells
probCIV = float(civ.count(1))/totalCells
probOVI = float(ovi.count(1))/totalCells

fo = open('probs.out', 'w')
s = 'Given \t Percent Probabiltiy of Finding\n'
fo.write(s)
s = '      \t HI \t MgII \t CIV \t OVI \n'
fo.write(s)
s = '{0:s} \t {1:>.3} \t {2:>.3} \t {3:>.3} \t {4:>.3} \n'.format('HI', probHI*100, probAgivenB(hi,mgii)*100, probAgivenB(hi,civ)*100, probAgivenB(hi,ovi)*100)
fo.write(s)
s = '{0:s} \t {1:>.3} \t {2:>.3} \t {3:>.3} \t {4:>.3} \n'.format('MgII', probAgivenB(mgii,hi)*100, probMgII*100, probAgivenB(mgii,civ)*100, probAgivenB(mgii,ovi)*100)
fo.write(s)
s = '{0:s} \t {1:>.3} \t {2:>.3} \t {3:>.3} \t {4:>.3} \n'.format('CIV',  probAgivenB(civ,hi)*100, probAgivenB(civ,mgii)*100, probCIV*100, probAgivenB(civ,ovi)*100)
fo.write(s)
s = '{0:s} \t {1:>.3} \t {2:>.3} \t {3:>.3} \t {4:>.3} \n'.format('OVI', probAgivenB(ovi,hi)*100, probAgivenB(ovi,mgii)*100, probAgivenB(ovi,civ)*100, probOVI*100)
fo.write(s)


fo.write('\n\n\n')

##########
# Calculate the probabilities for finding each ion if HI is NOT found
##########

# Loop through the files and pull out all ions that have absorption in any ion
# but not HI

f = open('dat.files')
'''
mgii_nH, mgii_t, mgii_z, mgii_l = [], [], [], []
civ_nH, civ_t, civ_z, civ_l = [], [], [], []
ovi_nH, ovi_t, ovi_z, ovi_l = [], [], [], []
mgiiCount = 0.
civCount = 0.
oviCount = 0.
for fN in f:
    ftmp = open(fN.strip())
    ftmp.readline()
    for line in ftmp:
        l = line.split()
        hiFlag = int(l[8])
        mgiiFlag = int(l[9])
        civFlag = int(l[10])
        oviFlag = int(l[11])

        
        if mgiiFlag==1 and hiFlag==0:
            mgii_nH.append( pow(10,float(l[4])))
            mgii_t.append(pow(10, float(l[5])))
            mgii_z.append(float(l[6]))
            mgii_l.append(float(l[7]))
            mgiiCount += 1
    
        if civFlag==1 and hiFlag==0:
            civ_nH.append(pow(10, float(l[4])))
            civ_t.append(pow(10, float(l[5])))
            civ_z.append(float(l[6]))
            civ_l.append(float(l[7]))
            civCount += 1

        if oviFlag==1 and hiFlag==0:
            ovi_nH.append(pow(10, float(l[4])))
            ovi_t.append(pow(10, float(l[5])))
            ovi_z.append(float(l[6]))
            ovi_l.append(float(l[7]))
            oviCount += 1

    ftmp.close()
f.close()

# Fraction of mgII cells that do not have hi to those that do
mgFrac = mgiiCount / mgii.count(1)
civFrac = civCount / civ.count(1)
oviFrac = oviCount / ovi.count(1)
'''

ion1, ion2, ion3 = findCellProps(f, 'HI')
mgii_nH = ion1[0]
mgii_t  = ion1[1]
mgii_z  = ion1[2]
mgii_l  = ion1[3]

civ_nH = ion2[0]
civ_t  = ion2[1]
civ_z  = ion2[2]
civ_l  = ion2[3]

ovi_nH = ion3[0]
ovi_t  = ion3[1]
ovi_z  = ion3[2]
ovi_l  = ion3[3]


print mgii_nH
print mgii_t
print mgii_z
print mgii_l

s = 'Cells that do not have HI absorption\n\n'
fo.write(s)
printCellProperties(fo, 'MgII', mgii_nH, mgii_t, mgii_z, mgii_l)
printCellProperties(fo, 'CIV', civ_nH, civ_t, civ_z, civ_l)
printCellProperties(fo, 'OVI', ovi_nH, ovi_t, ovi_z, ovi_l)

fo.write('\n\n\n')


# Repeat for cells that do not have MgII absorption










fo.close()
