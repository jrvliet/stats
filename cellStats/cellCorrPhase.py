

""" Plots the phase of gas cells and their correlationdd
    For example, plots the phase of cells that have CIV
    absorption given MgII, along with the cells that do 
    not have CIV given MgII
    """


import numpy as np
import matplotlib.pyplot as plt

def plotPhase(n, t, bins, ax, label):
    H, xedges, yedges = np.histogram2d(n, t, bins=bins)
    H = np.ma.masked_where(H==0, H)
    H = np.log10(H)
    H = np.rot90(H)
    H = np.flipud(H)
    mesh = ax.pcolormesh(xedges, yedges, H)
    ax.set_xlim([-8, 1])
    ax.set_ylim([2, 8])
    ax.text(-1,7,label)

    cbar = plt.colorbar(mesh, ax=ax, use_gridspec=True)
    cbar.ax.set_label('$\log$ (Counts) ')


numbins = 50

ions = ['HI', 'MgII', 'CIV', 'OVI']

c_given_mgii_nH = []
c_given_mgii_t = []
no_c_given_mgii_nH = []
no_c_given_mgii_t = []

mgii_given_c_nH = []
mgii_given_c_t = []
no_mgii_given_c_nH = []
no_mgii_given_c_t = []




n, t, hi, mgii, civ, ovi = [], [], [], [], [], []
# Read in the master list
f = open('cellMaster.dat')
f.readline()
for line in f:
    l = line.split()
    n.append(float(l[6]))
    t.append(float(l[7]))

    hi.append(int(l[10]))
    mgii.append(int(l[11]))
    civ.append(int(l[12]))
    ovi.append(int(l[13]))

f.close()
print 'Box read'
flags = [hi, mgii, civ, ovi]

# Loop over ion pairs
for i in range(0,len(ions)):

    for j in range(i+1,len(ions)):
        ion1_given_ion2_n = []
        ion1_given_ion2_t = []
        no_ion1_given_ion2_n = []
        no_ion1_given_ion2_t = []
        ion2_given_ion1_n = []
        ion2_given_ion1_t = []
        no_ion2_given_ion1_n = []
        no_ion2_given_ion1_t = []

        ion1 = ions[i]
        ion2 = ions[j]
        flags1 = flags[i]
        flags2 = flags[j]            

        print 'Ion1 = {0:s}\tIon2 = {1:s}'.format(ion1, ion2)
        for k in range(0,len(n)):
        
            f1 = flags1[k]
            f2 = flags2[k]
            if f1==1:
                if f2==1:
                    ion2_given_ion1_n.append(n[k])
                    ion2_given_ion1_t.append(t[k])
                else:
                    no_ion2_given_ion1_n.append(n[k])
                    no_ion2_given_ion1_t.append(t[k])
                    
            if f2==1:
                if f1==1:
                    ion1_given_ion2_n.append(n[k])
                    ion1_given_ion2_t.append(t[k])
                else:
                    no_ion1_given_ion2_n.append(n[k])
                    no_ion1_given_ion2_t.append(t[k])
                    
        fig, ((p11, p12), (p21, p22)) = plt.subplots(2,2,figsize=(12,12))
        plot_list = [p11, p21, p12, p22]
        labels = ['(a)', '(c)', '(b)', '(d)']

        # Plot 
        plotPhase(ion2_given_ion1_n, ion2_given_ion1_t, numbins, p11, '(a)')
        plotPhase(no_ion2_given_ion1_n, no_ion2_given_ion1_t, numbins, p12, '(b)')
        plotPhase(ion1_given_ion2_n, ion1_given_ion2_t, numbins, p21, '(c)')
        plotPhase(no_ion1_given_ion2_n, no_ion1_given_ion2_t, numbins, p22, '(d)')

        # Set titles
        p11.set_title(ion2+' given '+ion1)
        p12.set_title('No '+ion2+' given '+ion1)
        p21.set_title(ion1+'MgII given '+ion2)
        p22.set_title('No '+ion1+' given '+ion2)

        # Set x labels
        p21.set_xlabel('$\log (n_{H})$ [cm$^{-3}$] ')
        p22.set_xlabel('$\log (n_{H})$ [cm$^{-3}$] ')

        # Set y labels
        p11.set_ylabel(' $\log$ (T) [K] ')
        p21.set_ylabel(' $\log$ (T) [K] ')


        plt.tight_layout()
        s = ion1+'_'+ion2+'corrPhase.png'
        plt.savefig(s, bbox_inches='tight')



