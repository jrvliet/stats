# 
# This creates a list of all cells probed by an LOS
#


galID_list = ['D9o2', 'D9q', 'D9m4a']
expn_list1 = ['0.900', '0.926', '0.950', '0.975', '0.990', '1.002']
expn_list2 = ['0.901', '0.925', '0.950', '0.976', '0.991', '1.001']
expn_list3 = ['0.900', '0.925', '0.950', '0.975', '0.990', '1.000']
expn_list = [expn_list1, expn_list2, expn_list3]


base = '/home/jacob/matrix/sebass_gals/dwarfs/'
outputBase = '/home/jacob/matrix/sebass_gals/dwarfs/cellSummaries/'

base = '/home/matrix3/jrvander/sebass_gals/dwarfs/'
outputBase = '/home/matrix3/jrvander/sebass_gals/dwarfs/cellSummaries/'

for j in range(0,len(galID_list)):

    galID = galID_list[j]
    expn = expn_list[j]
#    print galID 
    for a in expn:
#        print '\t', a
        direc = base + galID + '_outputs/a'+a+'/'

        gasfile = direc+galID+'_GZa'+a+'.txt'
    
        f = open(gasfile)
        totalcells = 0
        for line in f:
            totalcells += 1
        f.close()

        direc += 'cellIDs/'

        ids = []

        for i in range(1,1000):
    
            losnum = '{0:d}'.format(i).zfill(4)
            filename = 'los'+losnum+'.cellID.dat'

            f = open(direc+filename)
            f.readline()
            for line in f:
                ids.append(int(line))

            f.close()

        uniqueIDs = set(ids)
        
        f = open(outputBase+galID+'.a'+a+'.probedLOS.dat', 'w')
        for i in uniqueIDs:
            f.write('{0:d}\n'.format(i))
        f.close()

       
        print 'For {0:s}, {1:s}: Probed {2:d} cells out of {3:d} ({4:%})'.format(galID, a, len(uniqueIDs), totalcells, len(uniqueIDs)/float(totalcells))

