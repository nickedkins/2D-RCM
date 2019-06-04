# Read and plot 2D RCM output data

# Just plot the ouput for the current run of PRRTM
import numpy as np
import matplotlib as matplotlib
import matplotlib.pyplot as plt
from pylab import *
from scipy import interpolate
from os import listdir
# import pandas as pd

directories = [
'_Current Output/'
]

# directories = [
'/Users/nickedkins/Dropbox/GitHub Repositories/Home/2D-RCM/_Useful Data/TOA res/tp=0.1/'
# ]


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
    boxlatcols = np.zeros((nlayersm,ncols))
    convcols = np.zeros((nlayersm,ncols))
 
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

    for col in range(ncols):    
        for i in range(nlayersm):
            boxlatcols[i,col] = f.readline()

    for col in range(ncols):    
        for i in range(nlayersm):
            convcols[i,col] = f.readline()

    sol_inc = 1362.0/4.0
    # abs_h2o = sum(abspncols[:nlayersm-1,:])*sol_inc/ncols / factor
    # abs_o3 = sum(A_oz_lcols[2:nlayersm,:])*sol_inc/ncols / factor
    # abs_surf = np.mean(abs_surf_lhcols[0,:]) / factor


    return tzmcols,pzmcols,wklm1cols,totuflumcols,htrmcols,altzmcols,pavelmcols,htro3cols,totdflumcols,wklm2cols,A_oz_lcols,abspncols,abs_surf_lhcols,tboundmcols,tavelmcols,nlayersm,ncols,boxlatcols,htrh2ocols,wklm3cols,convcols

# nlayersms=[31,100]
# ncols=5

plot_all_vert_profiles = 1

i1 = 0

tzm_store = np.zeros( ( 200, len(directories), 100 ) )
pico2_store = np.zeros( (len(directories), 100 ) )

tzm_master = []
pzm_master = []
boxlatcols_master = []

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
        tzmcols,pzmcols,wklm1cols,totuflumcols,htrmcols,altzmcols,pavelmcols,htro3cols,totdflumcols,wklm2cols,A_oz_lcols,abspncols,abs_surf_lhcols,tboundmcols,tavelmcols,nlayersm,ncols,boxlatcols,htrh2ocols,wklm3cols,convcols = readfile(fn,counter)
        tzm_master.append(tzmcols)  
        pzm_master.append(pzmcols)
        boxlatcols_master.append(boxlatcols)

        htrmlwcols = htrmcols[1:,:] - htro3cols - htrh2ocols

        i3=0

        if (plot_all_vert_profiles == 1):

            for col in range(ncols):

                conv_trop_ind = int(convcols[0,col])
                if conv_trop_ind > nlayersm:
                    conv_trop_ind = nlayersm

                p_trop = pzmcols[conv_trop_ind,col]
                t_trop = tzmcols[conv_trop_ind,col]
                z_trop = altzmcols[conv_trop_ind,col]

                print p_trop, t_trop, z_trop/1000.

                plt.figure(i2+1)
                plt.subplot(341)
                plt.title('tzm')
                plt.semilogy(tzmcols[:,col],pzmcols[:,col],'-o',label=str(fn))
                # plt.semilogy(tzmcols[:,col],pzmcols[:,col],label=str(fn))
                # plt.plot(tzmcols[:,col],altzmcols[:,col],'-o',label=str(fn))
                plt.plot(tzmcols[conv_trop_ind,col],pzmcols[conv_trop_ind,col],'*',markersize=20)
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))
                # plt.legend()
                
                plt.figure(i2+1)
                plt.subplot(342)
                plt.title('htrm')
                plt.semilogy(htrmcols[:,col],pzmcols[:,col],'-o')
                plt.xlim(-5,5)
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))
                plt.axvline(-0.03,ls='--')
                plt.axvline(0.03,ls='--')
                
                plt.figure(i2+1)
                plt.subplot(343)
                plt.title('totuflum')
                plt.semilogy(totuflumcols[:,col],pzmcols[:,col],'-o')
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))
                
                plt.figure(i2+1)
                plt.subplot(344)
                plt.title('totdflum')
                plt.semilogy(totdflumcols[:,col],pzmcols[:,col],'-o')
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))

                plt.figure(i2+1)
                plt.subplot(345)
                plt.title('abs_o3')
                plt.semilogy(A_oz_lcols[:,col]*1362./4.,pzmcols[1:,col],'-o')
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))      

                plt.figure(i2+1)
                plt.subplot(346)
                plt.title('abspn')
                plt.semilogy(abspncols[:,col]*1362./4.,pzmcols[1:,col],'-o')
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))

                plt.figure(i2+1)
                plt.subplot(347)
                plt.title('tavelm')
                plt.semilogy(tavelmcols[:,col],pzmcols[1:,col],'-o',label=str(fn))
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))
                # plt.legend()

                plt.figure(i2+1)
                plt.subplot(348)
                plt.title('htro3')
                plt.semilogy(htro3cols[:,col],pzmcols[1:,col],'-o',label=str(fn))
                plt.xlim(-5,5)
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))

                plt.figure(i2+1)
                plt.subplot(349)
                plt.title('htrh2o')
                plt.semilogy(htrh2ocols[:,col],pzmcols[1:,col],'-o',label=str(fn))
                plt.xlim(-5,5)
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))

                plt.figure(i2+1)
                plt.subplot(3,4,10)
                plt.title('wklm1 (h2o)')
                plt.loglog(wklm1cols[:,col],pzmcols[1:,col],'-o',label=str(fn))
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))

                plt.figure(i2+1)
                plt.subplot(3,4,11)
                plt.title('htrmlw')
                plt.semilogy(htrmlwcols[:,col],pzmcols[1:,col],'-o')
                plt.xlim(-5,5)
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))
                plt.axvline(-0.03,ls='--')
                plt.axvline(0.03,ls='--')

                plt.figure(i2+1)
                plt.subplot(3,4,12)
                plt.title('wklm3 (o3)')
                plt.semilogy(wklm3cols[:,col],pzmcols[1:,col],'-o',label=str(fn))
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))

                i3+=1

        # i2 += 1

    i1 += 1


# master indices: master[file][layer][column]
# tzm_master = np.array(tzm_master)
# pzm_master = np.array(pzm_master)
# boxlatcols_master = np.array(boxlatcols_master)

# filenames = np.array(filenames[0])

# x = np.linspace(-90,90,ncols)



# for i in range( shape(tzm_master)[0] ):

#     lats = boxlatcols_master[i,0,:]
#     pzms = pzm_master[i,:,0]
#     tzms = tzm_master[i,:,:]
#     tsgm_weighted = sum(tzms[0,:] * np.cos(np.deg2rad(lats)) / sum(np.cos(np.deg2rad(lats)) ))

#     plt.figure(1)
#     plt.subplot(121+i)
#     plt.contourf(lats,pzms,tzms,20)
#     plt.gca().set_yscale('log')
#     plt.ylim(1000,10)
#     plt.xlabel('Latitude')
#     plt.ylabel('Altitude')
#     plt.colorbar()

#     plt.figure(2)
#     plt.plot(lats,tzms[0,:])
    

# plt.subplot(223)
# plt.gca().set_yscale('log')
# plt.ylim(1000,10)
# plt.contourf(x,y,tzm_master[1,:,:]-tzm_master[0,:,:],20)
# plt.xlabel('Latitude')
# plt.ylabel('Altitude')
# plt.colorbar()

print 'Done'

############################################################
plt.tight_layout()
show()