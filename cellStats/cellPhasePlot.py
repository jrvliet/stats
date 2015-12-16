

import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt

numbins = 500

ion_list = ['HI', 'MgII', 'CIV', 'OVI']

# Read in cell data
nH   = []
t    = []
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
        nHTmp   = float(l2[4])
        tTmp    = float(l2[5])

        hi.append(hiTmp)
        mgii.append(mgiiTmp)
        civ.append(civTmp)
        ovi.append(oviTmp)
        nH.append(nHTmp)
        t.append(tTmp)

    ftmp.close()

f.close()

print np.min(hi)
print np.max(hi)
print np.min(mgii)
print np.max(mgii)
print np.min(civ)
print np.max(civ)
print np.min(ovi)
print np.max(ovi)

absorption = [hi, mgii, civ, ovi]
histos = []
xed = []
yed = []

numCells = len(hi)

# Bin up the data
for i in range(0,len(ion_list)):

    ion1 = ion_list[i]
    
    for j in range(0,len(ion_list)):

        ion2 = ion_list[j]

        density = []
        temperature = []
        
        # Collect the cells that have both ion1
        # and ion2 absorption
        print ion1, ion2
        for k in range(0,numCells):
            
            if absorption[i][k]==1 and  absorption[j][k]==1:
                density.append(nH[k])
                temperature.append(t[k])

        print len(density), len(temperature)
        
        H, xedges, yedges = np.histogram2d( density, temperature, bins=numbins)
        print np.mean(H)
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        Hmasked = np.ma.masked_where(H==0,H)
        Hmasked = np.log10(Hmasked)
        Hmasked = np.rot90(Hmasked)
        Hmasked = np.flipud(Hmasked)

        xed.append(xedges)
        yed.append(yedges)
        histos.append(Hmasked)


# Plot the data

fig,((p11,p12,p13,p14), 
     (p21,p22,p23,p24), 
     (p31,p32,p33,p34), 
     (p41,p42,p43,p44)) = plt.subplots(4,4,figsize=(12,12))

plot_list = [p11, p21, p31, p41, 
             p12, p22, p32, p42, 
             p13, p23, p33, p43, 
             p14, p24, p34, p44]


labels = ['(a)', '(e)', '(i)', '(m)', 
          '(b)', '(f)', '(j)', '(n)', 
          '(c)', '(g)', '(k)', '(o)',
          '(d)', '(h)', '(l)', '(p)'] 

for i in range(0,len(plot_list)):

    ax = plot_list[i]
    H = histos[i]
    xedges = xed[i]
    yedges = yed[i]
    mesh = ax.pcolormesh(xedges, yedges, H)
    ax.set_xlim([-8, 1])
    ax.set_ylim([2,8])
    ax.text(-1, 7, labels[i])
    
    if i==0:
        ax.set_title('HI', size=12)
        ax.set_ylabel('HI \n $\log$ (T) [K] ', multialignment='center')
    elif i==1:
        ax.set_ylabel('MgII \n $\log$ (T) [K] ', multialignment='center')
    elif i==2:
        ax.set_ylabel('CIV \n $\log$ (T) [K] ', multialignment='center')
    elif i==3:
        ax.set_ylabel('OVI \n $\log$ (T) [K] ', multialignment='center')
        ax.set_xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
    elif i==4:
        ax.set_title('MgII', size=12)
    elif i==7:
        ax.set_xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
    elif i==8:
        ax.set_title('CIV', size=12)
    elif i==11:
        ax.set_xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
    elif i==12:
        ax.set_title('OVI', size=12)
    elif i==15:
        ax.set_xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
    
    cbar = plt.colorbar(mesh, ax=ax, use_gridspec=True)
    cbar.ax.set_label('$\log$ (Counts)')


dropx = [p11, p12, p13, p14, p21, p22, p23, p24, p31, p32, p33, p34]
plt.setp([a.get_xticklabels() for a in dropx],visible=False)
       
dropy = [p12, p22, p32, p13, p23, p33, p42, p43, p14, p24, p34, p44]
plt.setp([a.get_yticklabels() for a in dropy],visible=False)

plt.tight_layout()

s = 'cellPhase.png'
plt.savefig(s, bbox_inches='tight')

