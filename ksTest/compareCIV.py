
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt

def shiftarr_leftedge(x,y):
    # From EWuncCDF.py
    # if first values start at 0,0
    
    # New arrays with twice the length
    newx = [None]*len(x)*2
    newy = [None]*len(y)*2

    # Start new array with first value
    #newx[0] = x[0]
    #newy[0] = y[0]
    
    k = 0

    # Continue with the rest of the values
    for j in range(len(y)):
        newx[k] = x[j]
        newy[k] = y[j]
        
        k += 1
        
        if j < (len(y)-1):
            newx[k] = x[j+1]
            newy[k] = y[j]
            k += 1

    newx[k] = x[-1]+5
    newy[k] = y[-1]
            
    return newx,newy


def mkcdf(x):
    # Make a cumulative distribution function out of an array of data
    # This array does not have to be sorted, this function sorts it

    # New arrays for sorted data and CDF value
    newx = [None]*(len(x)+1)
    newy = [None]*(len(x)+1)

    # Sort x array and place it into newx
    newx = np.sort(x)

    # Put a 0 at the start of the array
    newx = np.insert(newx,0,0)

    # Figure out the values of the cdf
    dy = 1.0/len(x)

    newy[0] = 0

    i = 1
    while i <= len(x):
        newy[i] = newy[i-1] + dy
        i += 1

    return newx,newy


# Read in COS-Dwarfs data
cosDataFileLoc = '/home/jacob/matrix/data/obs_data/cos_dwarfs/'
cosDataFile = 'cosdwarfs_detections.dat'
cosD = []
cosEW = []

f = open(cosDataFileLoc + cosDataFile)
for i in range(0,3):
    f.readline()
for line in f:

    l = line.split()
    r = float(l[7])
    rvir = float(l[8])
    cosD.append(r/rvir)
    cosEW.append(float(l[13])/1000.)

f.close()


# Read in Simulation Data
allFileLoc = '/home/jacob/matrix/sebass_gals/dwarfs/ALLfiles/masters/'

# Start with dwSN (D9o2)
allFile = 'D9o2.CIV.ALL.sysabs.large.master'
o2EW, o2D = np.loadtxt(allFileLoc+allFile, skiprows=1, usecols=(5, 22), unpack=True)

allFile = 'D9q.CIV.ALL.sysabs.large.master'
qEW, qD = np.loadtxt(allFileLoc+allFile, skiprows=1, usecols=(5, 22), unpack=True)

allFile = 'D9m4a.CIV.ALL.sysabs.large.master'
m4aEW, m4aD = np.loadtxt(allFileLoc+allFile, skiprows=1, usecols=(5, 22), unpack=True)


print o2EW[0], qEW[0], m4aEW[0]

# Perform the ks test
o2_stat,  o2_p  = stats.ks_2samp(cosEW, o2EW)
q_stat,   q_p   = stats.ks_2samp(cosEW, qEW)
m4a_stat, m4a_p = stats.ks_2samp(cosEW, m4aEW)


# Print the results
print 'dwSN: \n\t KS stat: {0:.5e}    \n\t p-value: {1:0.5e} \n\t percent: {2:.5e}\n'.format(o2_stat,  o2_p,  1.0-o2_p)
print 'dwALL_1: \n\t KS stat: {0:.5e} \n\t p-value: {1:0.5e} \n\t percent: {2:.5e}\n'.format(q_stat,   q_p,   1.0-q_p)
print 'dwALL_8: \n\t KS stat: {0:.5e} \n\t p-value: {1:0.5e} \n\t percent: {2:.5e}\n'.format(m4a_stat, m4a_p, 1.0-m4a_p)


print o2_stat-q_stat
print o2_stat-m4a_stat


# Plotting shenanigans

cosEWx, cosEWy = mkcdf(cosEW)
o2EWx,  o2EWy  = mkcdf(o2EW)
qEWx,   qEWy   = mkcdf(qEW)
m4aEWx, m4aEWy = mkcdf(m4aEW)

cosEWx, cosEWy = shiftarr_leftedge(cosEWx,cosEWy)
o2EWx,  o2EWy  = shiftarr_leftedge(o2EWx,o2EWy)
qEWx,   qEWy   = shiftarr_leftedge(qEWx,qEWy)
m4aEWx, m4aEWy = shiftarr_leftedge(m4aEWx,m4aEWy)

plt.plot(cosEWx, cosEWy, label='COS')
plt.plot(o2EWx,  o2EWy,  label='dwSN')
plt.plot(qEWx,   qEWy,   label='dwALL1')
plt.plot(m4aEWx, m4aEWy, label='dwALL8')

plt.legend(loc='lower right', frameon=False)
plt.xlim([0,1])
plt.savefig('play.pdf')



# Print results to file
fo = open('compareCIVResults.out', 'w')
fo.write( 'dwSN: \n\t KS stat: {0:.5e}    \n\t p-value: {1:0.5e} \n\t percent: {2:.5e}\n'.format(o2_stat,  o2_p,  1.0-o2_p) )
fo.write( 'dwALL_1: \n\t KS stat: {0:.5e} \n\t p-value: {1:0.5e} \n\t percent: {2:.5e}\n'.format(q_stat,   q_p,   1.0-q_p) )
fo.write( 'dwALL_8: \n\t KS stat: {0:.5e} \n\t p-value: {1:0.5e} \n\t percent: {2:.5e}\n'.format(m4a_stat, m4a_p, 1.0-m4a_p) )
fo.close()
