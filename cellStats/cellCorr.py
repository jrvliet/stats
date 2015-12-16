
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt


#civ = np.empty(1)
#ovi = np.empty(1)

#print civ

f = open('dat.files')

hi = []
mgii = []
civ = []
ovi = []

for line in f:

    print line
    l = line.strip()
    
    ftmp = open(l)
    ftmp.readline()
    for line2 in ftmp:
        
        l2 = line2.split()
        hiTmp = float(l2[5])
        mgiiTmp = float(l2[6])
        civTmp = float(l2[7])
        oviTmp = float(l2[8])
        
        hi.append(hiTmp)
        mgii.append(mgiiTmp)

        civ.append(civTmp)
        ovi.append(oviTmp)

    ftmp.close()
    print len(civ)
f.close()


print '\n\n File Reads Complete \n'

hiRelation = np.zeros(4)
mgiiRelation = np.zeros(4)
civRelation = np.zeros(4)
oviRelation = np.zeros(4)
foundCIV = 0
foundOVI = 0
foundBoth = 0
foundNone = 0
total = 0.
for i in range(0,len(civ)):

    total += 1
    h = hi[i]
    m = mgii[i]
    c = civ[i]
    o = ovi[i]
    
    if h==1:
        if m==1:
            hiRelation[1] += 1
        if c==1:
            hiRelation[2] += 1
        if o==1:
            hiRelation[3] += 1
        if m==0 and c==0 and o==0:
            hiRelation[0] += 1
    
    if m==1:
        if h==1:
            mgiiRelation[0] += 1
        if c==1:
            mgiiRelation[2] += 1
        if o==1:
            mgiiRelation[3] += 1
        if h==0 and c==0 and o==0:
            mgiiRelation[1] += 1
    
    if c==1:
        if m==1:
            civRelation[1] += 1
        if h==1:
            civRelation[0] += 1
        if o==1:
            civRelation[3] += 1
        if m==0 and c==0 and o==0:
            civRelation[2] += 1
    
    if o==1:
        if m==1:
            oviRelation[1] += 1
        if c==1:
            oviRelation[2] += 1
        if h==1:
            oviRelation[0] += 1
        if m==0 and c==0 and h==0:
            oviRelation[3] += 1
    

print '\n\n Loop Finished \n'

#founds = [foundNone, foundCIV, foundOVI, foundBoth]
bigList = [hiRelation, mgiiRelation, civRelation, oviRelation]
founds = []
for i in range(0,4):
    for j in range(0,4):
        val = bigList[i][j]
        if val==0:
            val=1
        founds.append(val)
print founds

#founds[3] = 1.
founds = np.log10(founds)
percents = []
#for found in founds:
#    per = float(found) / total
#    percents.append( '{0:.3%}'.format(per) )
## Make a bar plot of the results
N = 16
ind = np.arange(N)
width = 0.15

fig, ax = plt.subplots(1,1,figsize=(10.2,10.8))
rects1 = ax.bar(ind, founds, width )

ax.set_ylabel('Number of cells')
ax.set_xticks(ind)
ax.set_xticklabels( ('Only HI', 'HI and MgII', 'HI and CIV', 'HI and OVI', 'MgII and HI', 'Only MgII', 'MgII and CIV', 'MgII and OVI', 'CIV and HI', 'CIV and MgII', 'Only CIV', 'CIV and OVI', 'OVI and HI', 'OVI and MgII', 'OVI and CIV', 'Only OVI') )

plt.setp(ax.xaxis.get_majorticklabels(), rotation=90, size=10)
plt.tick_params(axis='x', which='both', bottom='off', top='off')
plt.ylim([3,6])
def autolabel(rects, values):
    for i in range(0,len(rects)):
        rect = rects[i]
        val = values[i]
        print val
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, val[i],
                ha='center', va='bottom')

#autolabel(rects1, percents)

plt.savefig('cells.png')

#p = st.spearmanr(civ, ovi)
#print p
