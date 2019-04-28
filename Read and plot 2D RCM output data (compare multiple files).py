# Read and plot 2D RCM output data

# Just plot the ouput for the current run of PRRTM
import numpy as np
import matplotlib as matplotlib
import matplotlib.pyplot as plt
from pylab import *
from scipy import interpolate
from os import listdir
# import pandas as pd

# directories = [
# '_Current Output/'
# ]

directories = [
'_Useful Data/lapse rate varied/'
]

linestyles = ['-','--','--']

colors = ['b','r','g','orange','purple','yellow','pink']

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
            tzmcols[i,col] = f.readline()

    for col in range(ncols):
        for i in range(nlayersm+1):
            altzmcols[i,col] = f.readline()

    for col in range(ncols):
        for i in range(nlayersm+1):
            totuflumcols[i,col] = f.readline()

    for col in range(ncols):
        for i in range(nlayersm+1):
            totdflumcols[i,col] = f.readline()

    for col in range(ncols):
        for i in range(nlayersm+1):
            htrmcols[i,col] = f.readline()

    for col in range(ncols):
        for i in range(nlayersm):
            htrh2ocols[i,col] = f.readline()

    for col in range(ncols):
        for i in range(nlayersm):
            htro3cols[i,col] = f.readline()

    for col in range(ncols):
        for i in range(nlayersm):
            wklm1cols[i,col] = f.readline()

    for col in range(ncols):
        for i in range(nlayersm):
            wklm2cols[i,col] = f.readline()

    for col in range(ncols):
        for i in range(nlayersm):
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

    sol_inc = 1362.0/4.0
    # abs_h2o = sum(abspncols[:nlayersm-1,:])*sol_inc/ncols / factor
    # abs_o3 = sum(A_oz_lcols[2:nlayersm,:])*sol_inc/ncols / factor
    # abs_surf = np.mean(abs_surf_lhcols[0,:]) / factor


    return tzmcols,pzmcols,wklm1cols,totuflumcols,htrmcols,altzmcols,pavelmcols,htro3cols,totdflumcols,wklm2cols,A_oz_lcols,abspncols,abs_surf_lhcols,tboundmcols,tavelmcols,nlayersm,ncols

# nlayersms=[31,100]
# ncols=5

plot_all_vert_profiles = 1

i1 = 0

tzm_store = np.zeros( ( 200, len(directories), 100 ) )
pico2_store = np.zeros( (len(directories), 100 ) )

tzm_master = []
pzm_master = []

filenames = []

for directory in directories:

    ls = linestyles[i1]
    color = colors[i1]

    a = sorted(listdir(directory))

    filenames.append(a)

    counter=0

    i2 = 0

    for fn in a:
        if (fn == '.DS_Store' or fn == 'new benchmark'):
            continue
        tzmcols,pzmcols,wklm1cols,totuflumcols,htrmcols,altzmcols,pavelmcols,htro3cols,totdflumcols,wklm2cols,A_oz_lcols,abspncols,abs_surf_lhcols,tboundmcols,tavelmcols,nlayersm,ncols = readfile(fn,counter)
        tzm_master.append(tzmcols)  
        pzm_master.append(pzmcols)

        i3=0

        if (plot_all_vert_profiles == 1):

            for col in range(ncols):

                plt.figure(1)
                plt.subplot(331)
                plt.title('tzm')
                plt.semilogy(tzmcols[:,col],pzmcols[:,col],label=str(fn),c=colors[i3],ls=linestyles[i2])
                # plt.semilogy(tzmcols[:,col],pzmcols[:,col],label=str(fn))
                # plt.plot(tzmcols[:,col],altzmcols[:,col],'-o',label=str(fn))
                plt.plot(tboundmcols[0,col],pzmcols[0,col],'*',c=colors[i3],markersize=20)
                plt.ylim(pzmcols[0,col]*1.1,1)
                # plt.legend()
                
                plt.subplot(332)
                plt.title('htrm')
                plt.semilogy(htrmcols[:,col],pzmcols[:,col],ls=linestyles[i2],c=colors[i3])
                plt.ylim(pzmcols[0,col]*1.1,1)
                plt.axvline(-0.03,ls='--')
                plt.axvline(0.03,ls='--')
                
                plt.subplot(333)
                plt.title('totuflum')
                plt.semilogy(totuflumcols[:,col],pzmcols[:,col],ls=linestyles[i2],c=colors[i3])
                plt.ylim(pzmcols[0,col]*1.1,1)
                
                plt.subplot(334)
                plt.title('totdfllum')
                plt.semilogy(totdflumcols[:,col],pzmcols[:,col],ls=linestyles[i2],c=colors[i3])
                plt.ylim(pzmcols[0,col]*1.1,1)

                plt.subplot(335)
                plt.title('abs_o3')
                plt.semilogy(A_oz_lcols[:,col]*1362./4.,pzmcols[1:,col],ls=linestyles[i2],c=colors[i3])
                plt.ylim(pzmcols[0,col]*1.1,1)      

                plt.subplot(336)
                plt.title('abspn')
                plt.semilogy(abspncols[:,col]*1362./4.,pzmcols[1:,col],ls=linestyles[i2],c=colors[i3])
                plt.ylim(pzmcols[0,col]*1.1,1)

                plt.subplot(337)
                plt.title('tavelm')
                plt.semilogy(tavelmcols[:,col],pzmcols[1:,col],label=str(fn),c=colors[i3],ls=linestyles[i2])
                plt.ylim(pzmcols[0,col]*1.1,1)
                # plt.legend()

                i3+=1

        i2 += 1

    i1 += 1

# master indices: master[file][layer][column]

tzm_master = np.array(tzm_master)
pzm_master = np.array(pzm_master)

filenames = np.array(filenames[0])

# for file in range(len(a)):
#     if(filenames[file]=='baseline'):
#         plt.plot(tzm_master[file,0,:],label=filenames[file],linestyle='--')    
#     else:
#         plt.plot(tzm_master[file,0,:],label=filenames[file])
# plt.legend()

# print np.mean( tzm_master[ :,0,: ], axis=1 )

lapses = np.linspace(0,10,6)
# plt.plot(lapses, np.mean( tzm_master[ :,0,:], axis=1 ) - np.mean( tzm_master[ 0,0,: ], axis=0 ),'-o' )
# plt.plot(lapses, np.mean( tzm_master[ :,0,:], axis=1 ),'-o' )
# plt.plot( lapses,tzm_master[ :,0,0],'-o')
# plt.xlabel('Polar box lapse rate (K/km)')
# plt.ylabel(r'$\Delta T$ compared to $\Gamma=0$')
    


print 'Done'

############################################################
plt.tight_layout()
show()