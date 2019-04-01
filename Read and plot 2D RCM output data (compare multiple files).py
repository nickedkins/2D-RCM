# Read and plot 2D RCM output data

# Just plot the ouput for the current run of PRRTM
import numpy as np
import matplotlib as matplotlib
import matplotlib.pyplot as plt
from pylab import *
from scipy import interpolate
from os import listdir
# from pandas import *

directories = [
'_Current Output/'
]

# directories = [
# '_Useful Data/vary sa/',
# # '_Useful Data/vary sa tboundm tweak/',
# ]


linestyles = ['-','--','-o']

colors = ['b','r','g']

def init_plotting():
    plt.rcParams['figure.figsize'] = (10,10)
    plt.rcParams['font.size'] = 15
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['axes.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['axes.titlesize'] = 1.2*plt.rcParams['font.size']
    plt.rcParams['legend.fontsize'] = plt.rcParams['font.size']
    plt.rcParams['xtick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['ytick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['savefig.dpi'] = 2*plt.rcParams['savefig.dpi']
    plt.rcParams['xtick.major.size'] = 3    
    plt.rcParams['xtick.minor.size'] = 3
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.minor.width'] = 1
    plt.rcParams['ytick.major.size'] = 3
    plt.rcParams['ytick.minor.size'] = 3
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['ytick.minor.width'] = 1
    plt.rcParams['legend.frameon'] = True
    plt.rcParams['legend.loc'] = 'best'
    plt.rcParams['axes.linewidth'] = 1

    plt.rcParams['lines.linewidth'] = 1.0   
    plt.rcParams['lines.markersize'] = 8 
    plt.rcParams['figure.facecolor'] = 'white'
    plt.rcParams['axes.facecolor'] = 'white'
    #plt.rcParams['axes.color_cycle'] = ['b', 'r', 'g','pink','orange','darkgreen','purple']

    plt.rcParams['grid.color'] = 'k'
    plt.rcParams['grid.linestyle'] = ':'
    plt.rcParams['grid.linewidth'] = 0.5

    #plt.gca().spines['right'].set_color('None')
    #plt.gca().spines['top'].set_color('None')
    plt.gca().xaxis.set_ticks_position('bottom')
    plt.gca().yaxis.set_ticks_position('left')
init_plotting()


def readfile(fn,counter):
    # print str(fn)

    output_file = directory + fn

    f = open(output_file,'r')
    

    params = f.readline().split()

    ncols=int(params[0])
    nlayersm=int(params[1])

    tzmcols = np.zeros((nlayersm+1,ncols))
    altzmcols = np.zeros((nlayersm+1,ncols))
    totuflumcols = np.zeros((nlayersm+1,ncols))
    totdflumcols = np.zeros((nlayersm+1,ncols))
    htrmcols = np.zeros((nlayersm+1,ncols))
    htrh2ocols = np.zeros((nlayersm,ncols))
    htro3cols = np.zeros((nlayersm,ncols))
    wklm1cols = np.zeros((nlayersm,ncols))
    wklm2cols = np.zeros((nlayersm,ncols))
    wklm3cols = np.zeros((nlayersm,ncols))
    wbrodlmcols = np.zeros((nlayersm,ncols))
    abspncols = np.zeros((nlayersm,ncols))
    A_oz_lcols = np.zeros((nlayersm,ncols))
    abs_surf_lhcols = np.zeros((nlayersm,ncols))
    altlaymcols = np.zeros((nlayersm,ncols))
    pzmcols = np.zeros((nlayersm+1,ncols))
    pavelmcols = np.zeros((nlayersm,ncols))
    tavelmcols = np.zeros((nlayersm,ncols))
    tboundmcols= np.zeros((nlayersm,ncols))
    R_gcols = np.zeros((nlayersm,ncols))
 


    for col in range(ncols):
        for i in range(nlayersm+1):
            # print col, i, f.readline()
            tzmcols[i,col] = f.readline()

    for col in range(ncols):
        for i in range(nlayersm+1):
            # print col, i, f.readline()
            altzmcols[i,col] = f.readline()

    for col in range(ncols):
        for i in range(nlayersm+1):
            # print col, i, f.readline()
            totuflumcols[i,col] = f.readline()

    for col in range(ncols):
        for i in range(nlayersm+1):
            # print col, i, f.readline()
            totdflumcols[i,col] = f.readline()

    for col in range(ncols):
        for i in range(nlayersm+1):
            # print col, i, f.readline()
            htrmcols[i,col] = f.readline()

    for col in range(ncols):
        for i in range(nlayersm):
            # print col, i, f.readline()
            htrh2ocols[i,col] = f.readline()

    for col in range(ncols):
        for i in range(nlayersm):
            # print col, i, f.readline()
            htro3cols[i,col] = f.readline()

    for col in range(ncols):
        for i in range(nlayersm):
            # print col, i, f.readline()
            wklm1cols[i,col] = f.readline()

    for col in range(ncols):
        for i in range(nlayersm):
            # print col, i, f.readline()
            wklm2cols[i,col] = f.readline()

    for col in range(ncols):
        for i in range(nlayersm):
            # print col, i, f.readline()
            wklm3cols[i,col] = f.readline()

    for col in range(ncols):
        for i in range(nlayersm):
            wbrodlmcols[i,col] = f.readline()

    for col in range(ncols):
        for i in range(nlayersm):
            abspncols[i,col] = f.readline()

    for col in range(ncols):
        for i in range(nlayersm):
            A_oz_lcols[i,col] = f.readline()

    for col in range(ncols):    
        for i in range(nlayersm):
            abs_surf_lhcols[i,col] = f.readline()

    for col in range(ncols):    
        for i in range(nlayersm):
            altlaymcols[i,col] = f.readline()

    for col in range(ncols):    
        for i in range(nlayersm+1):
            pzmcols[i,col] = f.readline()            

    for col in range(ncols):    
        for i in range(nlayersm):
            pavelmcols[i,col] = f.readline()

    for col in range(ncols):    
        for i in range(nlayersm):
            tavelmcols[i,col] = f.readline()

    for col in range(ncols):    
        for i in range(nlayersm):
            tboundmcols[i,col] = f.readline()

    for col in range(ncols):    
        for i in range(nlayersm):
            R_gcols[i,col] = f.readline()


    # print 'Total abs surf: ', abs_surf_lhcols[0,:], np.mean(abs_surf_lhcols[0,:])

    sol_inc = 1362.0/4.0


    # abs_h2o = sum(abspncols[:nlayersm-1,:])*sol_inc/ncols / factor
    # abs_o3 = sum(A_oz_lcols[2:nlayersm,:])*sol_inc/ncols / factor
    # abs_surf = np.mean(abs_surf_lhcols[0,:]) / factor


    return tzmcols,pzmcols,wklm1cols,totuflumcols,htrmcols,altzmcols,pavelmcols,htro3cols,totdflumcols,wklm2cols,A_oz_lcols,abspncols,abs_surf_lhcols,tboundmcols,tavelmcols

nlayersms=[31,31]
ncols=1



i1 = 0
i2 = 0

tzm_store = np.zeros( ( 200, len(directories), 100 ) )
pico2_store = np.zeros( (len(directories), 100 ) )

for directory in directories:

    nlayersm = nlayersms[i1]


    p = np.linspace(1000,1000.0/nlayersm,nlayersm)
    p1 = np.linspace(1000,1000.0/(nlayersm+1),nlayersm+1)

    ls = linestyles[i1]
    color = colors[i1]

    a = sorted(listdir(directory))

    print a

    counter=0

    for fn in a:
        if (fn == '.DS_Store' or fn == 'new benchmark'):
            continue
        tzmcols,pzmcols,wklm1cols,totuflumcols,htrmcols,altzmcols,pavelmcols,htro3cols,totdflumcols,wklm2cols,A_oz_lcols,abspncols,abs_surf_lhcols,tboundmcols,tavelmcols = readfile(fn,counter)
        pico2 = (pzmcols[0] - 800.0)/1000.0
        # pico2 = (pzmcols[0]/2.0)/1000.0

        
        plt.figure(1)
        plt.semilogy(tzmcols,pzmcols,ls=ls,label='{:3.0f} bar'.format(pico2[0]))
        plt.ylim(pzmcols[0]*1.1,1)
        plt.ylabel('Pressure (hPa)')
        plt.xlabel('Temperature (K)')
        plt.legend(title='$p_{I,N_2}$')

        plt.figure(1)
        plt.semilogy(totdflumcols-totuflumcols,pzmcols,ls=ls,label='{:3.0f} bar'.format(pico2[0]))
        plt.ylim(pzmcols[0]*1.1,1)
        plt.ylabel('Pressure (hPa)')
        plt.xlabel('Net Flux (Wm$^{-2}$)')
        plt.legend(title='$p_{I,N_2}$')

        plt.figure(1)
        plt.subplot(331)
        plt.title('tzm')
        # plt.semilogy(tzmcols,pzmcols,'-o',label=pico2)
        plt.plot(tzmcols,altzmcols,'-o',label=str(fn))
        # plt.ylim(pzmcols[0]*1.1,1)
        plt.legend()
        
        plt.subplot(332)
        plt.title('htrm')
        plt.semilogy(htrmcols,pzmcols,'-o',ls=ls)
        plt.ylim(pzmcols[0]*1.1,1)
        plt.axvline(-0.03,ls='--')
        plt.axvline(0.03,ls='--')
        
        plt.subplot(333)
        plt.title('totuflum')
        plt.semilogy(totuflumcols,pzmcols,'-o',ls=ls)
        plt.ylim(pzmcols[0]*1.1,1)
        
        plt.subplot(334)
        plt.title('totdfllum')
        plt.semilogy(totdflumcols,pzmcols,'-o',ls=ls)
        plt.ylim(pzmcols[0]*1.1,1)

        plt.subplot(335)
        plt.title('abs_o3')
        plt.semilogy(A_oz_lcols*1362./4.,pzmcols[1:],'-o',ls=ls)
        plt.ylim(pzmcols[0]*1.1,1)        

        plt.subplot(336)
        plt.title('abspn')
        plt.semilogy(abspncols*1362./4.,pzmcols[1:],'-o',ls=ls)
        plt.ylim(pzmcols[0]*1.1,1)

        plt.subplot(337)
        plt.title('tavelm')
        plt.semilogy(tavelmcols,pzmcols[1:],'-o',label=str(fn))
        plt.ylim(pzmcols[0]*1.1,1)
        plt.legend()


        # plt.figure(1)
        # # plt.subplot(337)
        # plt.title(r'Inversion develops with increased $\alpha_{surf}$')
        # plt.semilogy(tavelmcols,pzmcols[1:],'-o',label=str(fn))
        # plt.ylim(pzmcols[0]*1.1,1)
        # plt.xlabel('Temperature (K)')
        # plt.ylabel('Pressure (hPa)')
        # plt.legend()

        # print sum(abspncols*1362./4.), sum(A_oz_lcols*1362./4.), abs_surf_lhcols[1]

        # for i in range(nlayersm-1):
        #     print pzmcols[i][0],',',abspncols[i][0]*1362./4.,',',A_oz_lcols[i][0]*1362./4.
        
        # plt.figure(1)
        # plt.subplot(121)
        # plt.title('tzm')
        # plt.semilogy(tzmcols,pzmcols,'-',label=str(fn),ls=ls)
        # # plt.plot(pzmcols,altzmcols,'-o',label='%4.2f bar' % (pico2))
        # plt.xlabel('Temperature (K)')
        # plt.ylabel('Pressure (hPa)')
        # plt.ylim(pzmcols[0]*1.1,1)
        # plt.legend()
        # plt.grid(which='both')
        
        # plt.subplot(122)
        # plt.title('totuflum')
        # plt.semilogy(totdflumcols-totuflumcols,pzmcols,'-',label=str(fn),ls=ls)
        # # plt.plot(totuflumcols,altzmcols,'-o',label='%4.2f bar' % (pico2))
        # plt.xlabel('Net Downward Flux (Wm$^{-2}$)')
        # plt.ylabel('Pressure (hPa)')
        # plt.ylim(pzmcols[0]*1.1,1)
        # plt.legend()
        # plt.grid(which='both')
        

        # pico2_store[i1,i2] = pico2

        # tzm_store[:nlayersm,i1,i2] = tzmcols[:,0]

        # pico2_wklm = wklm2cols[1,0]/4.649531e+23
        
        # plt.figure(2)
        # plt.semilogx(pico2_wklm,tzmcols[0,0],'o',c=color)
        # plt.xlabel('CO$_2$ Inventory (bar)')
        # plt.ylabel('OLR (Wm$^{-2}$)')

        
        # plt.imshow(tzmcols,interpolation='none')

        # plt.plot(tzmcols[0,:],label='tzm[0]')
        # plt.plot(tavelmcols[1,:],label='tavelm[1]')
        # plt.legend()

        # plt.semilogy(tavelmcols[:,0],pzmcols[1:,0],'-o')
        # plt.ylim( max(pzmcols[:,2]), min(pzmcols[:,2]) )

        i2 += 1

    i1 += 1



# plt.figure(1)
# plt.semilogx(pico2_store[0,:],tzm_store[0,0,:],'-o',label='with PB')
# plt.semilogx(pico2_store[1,:],tzm_store[0,1,:],'-o',label='without PB')
# plt.xlabel('N$_2$ Inventory (bar)')
# plt.ylabel('Surface Temperature (K)')    
# plt.legend()

plt.tight_layout()
print 'Done'
show()