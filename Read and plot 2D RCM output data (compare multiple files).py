# Read and plot 2D RCM output data

# Just plot the ouput for the current run of PRRTM
import numpy as np
# import matplotlib as matplotlib
import matplotlib.pyplot as plt
from pylab import *
from scipy import interpolate
from os import listdir
import matplotlib.colors as colors
# import pandas as pd

plot_all_vert_profiles = 1
legends_on = 0
grids_on = 1

directories = [
'_Current Output/',
]

# set the colormap and centre the colorbar
class MidpointNormalize(colors.Normalize):
    """
    Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

    e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

def init_plotting():
    plt.rcParams['figure.figsize'] = (10,10)
    plt.rcParams['font.size'] = 15
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['axes.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['axes.titlesize'] = 1.2*plt.rcParams['font.size']
    plt.rcParams['legend.fontsize'] = plt.rcParams['font.size']
    plt.rcParams['xtick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['ytick.labelsize'] = plt.rcParams['font.size']
    # plt.rcParams['savefig.dpi'] = 2*plt.rcParams['savefig.dpi']
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

def latwghtavg(x,lats):
    # print "x:", x
    # print "lats:", lats
    # print "sinlats:", np.cos(np.deg2rad(lats))
    # print "sumlats:", sum(np.cos(np.deg2rad(lats)))
    x_avg = sum(x * np.cos(np.deg2rad(lats))) / sum(np.cos(np.deg2rad(lats)))
    return x_avg

def latwghtavg_2d(x,lats):
    nlays = len(x[:,0])
    x_avg = np.zeros(nlays)
    for i in range(nlays):
        x_avg[i] = sum(x[i,:] * np.cos(np.deg2rad(lats))) / sum(np.cos(np.deg2rad(lats)))
    return x_avg


linestyles = ['-','--','--','o']

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
    lapsecritcols = np.zeros((nlayersm,ncols))
    meridtransp = np.zeros((nlayersm,ncols))
    abs_h2o_cols = np.zeros((nlayersm,ncols))
    abs_o3_cols = np.zeros((nlayersm,ncols))
    abs_surf_cols = np.zeros((nlayersm,ncols))
    d_mid = np.zeros((nlayersm,ncols))
    d_trop = np.zeros((nlayersm,ncols))
    rel_hum_cols = np.zeros((nlayersm,ncols))

 
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

    for col in range(ncols):        
        for i in range(nlayersm):
            lapsecritcols[i,col] = f.readline()

    for col in range(ncols):        
        for i in range(nlayersm):
            meridtransp[i,col] = f.readline()

    for col in range(ncols):        
        for i in range(nlayersm):
            abs_h2o_cols[i,col] = f.readline()

    for col in range(ncols):        
        for i in range(nlayersm):
            abs_o3_cols[i,col] = f.readline()            

    for col in range(ncols):        
        for i in range(nlayersm):
            abs_surf_cols[i,col] = f.readline()

    for col in range(ncols):        
        for i in range(nlayersm):
            d_mid[i,col] = f.readline()

    for col in range(ncols):        
        for i in range(nlayersm):
            d_trop[i,col] = f.readline()            

    for col in range(ncols):        
        for i in range(nlayersm):
            rel_hum_cols[i,col] = f.readline()            


    sol_inc = 1362.0/4.0
    # abs_h2o = sum(abspncols[:nlayersm-1,:])*sol_inc/ncols / factor
    # abs_o3 = sum(A_oz_lcols[2:nlayersm,:])*sol_inc/ncols / factor
    # abs_surf = np.mean(abs_surf_lhcols[0,:]) / factor


    return tzmcols,pzmcols,wklm1cols,totuflumcols,htrmcols,altzmcols,pavelmcols,htro3cols,totdflumcols,wklm2cols,A_oz_lcols,\
    abspncols,abs_surf_lhcols,tboundmcols,tavelmcols,nlayersm,ncols,boxlatcols,htrh2ocols,wklm3cols,convcols,wbrodlmcols,\
    lapsecritcols,meridtransp,abs_h2o_cols,abs_o3_cols,abs_surf_cols,d_mid,d_trop,rel_hum_cols

i1 = 0

i_dir=0
for directory in directories:

    tzm_master = []
    pzm_master = []
    altzm_master = []
    boxlatcols_master = []
    lapsecritcols_master = []
    conv_trop_ind_master = []
    meridtransp_master = []
    abspncols_master = []
    A_oz_lcols_master = []
    abs_surf_lhcols_master = []
    totuflumcols_master = []
    abs_h2o_cols_master = []
    abs_o3_cols_master = []
    abs_surf_cols_master = []
    d_mid_master = []
    d_trop_master = []
    rel_hum_master = []
    wklm1_master = []
    wbrodlm_master = []

    filenames = []

    dir_label = directory.split('/')[-2]
    print(dir_label)

    ls = linestyles[i1]
    # color = colors[i1]

    a = sorted(listdir(directory))
    if('.DS_Store' in a):
        a.remove('.DS_Store')

    filenames.append(a)

    counter=0

    i2 = 0

    for fn in a:
        if (fn=='.DS_Store'):
            continue
        tzmcols,pzmcols,wklm1cols,totuflumcols,htrmcols,altzmcols,pavelmcols,htro3cols,totdflumcols,wklm2cols,A_oz_lcols,abspncols,\
        abs_surf_lhcols,tboundmcols,tavelmcols,nlayersm,ncols,boxlatcols,htrh2ocols,wklm3cols,convcols,wbrodlmcols,lapsecritcols,\
        meridtransp,abs_h2o_cols,abs_o3_cols,abs_surf_cols,d_mid,d_trop,rel_hum_cols = readfile(fn,counter)


        tzm_master.append(tzmcols)  
        pzm_master.append(pzmcols)
        altzm_master.append(altzmcols)  
        boxlatcols_master.append(boxlatcols)
        lapsecritcols_master.append(lapsecritcols)
        meridtransp_master.append(meridtransp)
        abspncols_master.append(abspncols)
        A_oz_lcols_master.append(A_oz_lcols)
        abs_surf_lhcols_master.append(abs_surf_lhcols)
        totuflumcols_master.append(totuflumcols)
        abs_h2o_cols_master.append(abs_h2o_cols)        
        abs_o3_cols_master.append(abs_o3_cols)        
        abs_surf_cols_master.append(abs_surf_cols)
        d_mid_master.append(d_mid)       
        d_trop_master.append(d_trop)
        rel_hum_master.append(rel_hum_cols)
        wklm1_master.append(wklm1cols)
        wbrodlm_master.append(wbrodlmcols)


        htrmlwcols = htrmcols[1:,:] - htro3cols - htrh2ocols

        conv_trop_ind_cols = np.zeros(ncols)

        for col in range(ncols):
            conv_trop_ind = int(convcols[0,col])
            if conv_trop_ind > nlayersm:    
                conv_trop_ind = nlayersm
            conv_trop_ind_cols[col] = int(conv_trop_ind)


        conv_trop_ind_cols = conv_trop_ind_cols.astype(int)
        conv_trop_ind_master.append(conv_trop_ind_cols)


        # plt.figure(1)

        # plt.subplot(221)
        # plt.plot(boxlatcols[1,:],-1.0*lapsecritcols[1,:],'-o',label=str(fn))
        # plt.xlabel('Latitude')
        # plt.ylabel('Lapse rate (K/km)')
        # plt.legend()

        # plt.subplot(222)
        # plt.plot(boxlatcols[1,:],tzmcols[0,:],'-o',label=str(fn))
        # plt.xlabel('Latitude')
        # plt.ylabel('Surface temperature (K)')
        # plt.legend()

        # for col in range(ncols):
        #     if(col==1):
        #         plt.subplot(223)
        #         plt.plot(boxlatcols[1,col],tzmcols[conv_trop_ind_cols[col],col],'-o',label=str(fn),color=colors[i2])
        #         plt.xlabel('Latitude')
        #         plt.ylabel('Tropopause temperature (K)')
        #         plt.legend()
        #     else:
        #         plt.subplot(223)
        #         plt.plot(boxlatcols[1,col],tzmcols[conv_trop_ind_cols[col],col],'-o',color=colors[i2])
        #         plt.xlabel('Latitude')
        #         plt.ylabel('Tropopause temperature (K)')
        #         plt.legend()

        # for col in range(ncols):
        #     if(col==1):
        #         plt.subplot(224)
        #         plt.plot(boxlatcols[1,col],altzmcols[conv_trop_ind_cols[col],col]/1000.,'-o',label=str(fn),color=colors[i2])
        #         plt.xlabel('Latitude')
        #         plt.ylabel('Tropopause altitude (km)')
        #         plt.legend()
        #     else:
        #         plt.subplot(224)
        #         plt.plot(boxlatcols[1,col],altzmcols[conv_trop_ind_cols[col],col]/1000.,'-o',color=colors[i2])
        #         plt.xlabel('Latitude')
        #         plt.ylabel('Tropopause altitude (km)')
        #         plt.legend()

        # plt.subplot(221)
        # plt.plot(boxlatcols[1,:],-1.0*lapsecritcols[1,:],'-o',label=str(fn))
        # plt.xlabel('Latitude')
        # plt.ylabel('Lapse rate (K/km)')
        # plt.legend()


        i3=0

        if (plot_all_vert_profiles == 1):

            for col in range(ncols):

                conv_trop_ind = int(convcols[0,col])
                if conv_trop_ind > nlayersm:    
                    conv_trop_ind = nlayersm

                # print('DW LW at tropopause: ', totdflumcols[conv_trop_ind,col])
                # print('abs SW above tropopause: ', sum(A_oz_lcols[conv_trop_ind:nlayersm,col])*1362./4.)

                p_trop = pzmcols[conv_trop_ind,col]
                t_trop = tzmcols[conv_trop_ind,col]
                z_trop = altzmcols[conv_trop_ind,col]

                print(p_trop,',', t_trop,',', z_trop/1000.,',', pzmcols[0,col],',', tzmcols[0,col])

                # plt.figure(2)
                # plt.plot(lapsecritcols[0,col],z_trop/1000.,'o')

                # for i in range(len(altzmcols[:,col])):
                #     print altzmcols[i,col]/1000., ',', pzmcols[i,col], ',', tzmcols[i,col]

                plt.figure(i1+1)
                plt.subplot(231)
                plt.title('tzm')
                plt.semilogy(tzmcols[:,col],pzmcols[:,col],label=str(fn))
                plt.plot(tzmcols[conv_trop_ind,col],pzmcols[conv_trop_ind,col],'*',markersize=20)
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))
                if(grids_on==2):
                    plt.gca().minorticks_on()
                if(grids_on==1):
                    plt.grid(which='both',axis='both')
                if(legends_on==1):
                    plt.legend()    
                
                plt.figure(i1+1)
                plt.subplot(232)
                plt.title('htrm')
                plt.semilogy(htrmcols[:,col],pzmcols[:,col],ls=linestyles[i1],label=str(fn))
                plt.xlim(-5,5)
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))
                plt.axvline(-0.1,ls='--')
                plt.axvline(0.1,ls='--')
                if(legends_on==1):
                    plt.legend()
                
                plt.figure(i1+1)
                plt.subplot(233)
                plt.title('totuflum')
                plt.semilogy(totuflumcols[:,col],pzmcols[:,col],ls=linestyles[i1],label='up')
                plt.semilogy(totdflumcols[:,col],pzmcols[:,col],ls=linestyles[i1],label='down')
                plt.semilogy(totdflumcols[:,col]-totuflumcols[:,col],pzmcols[:,col],ls=linestyles[i1],label='net')
                plt.axvline(0,ls='--')
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))
                if(legends_on==1):
                    plt.legend()

                plt.figure(i1+1)
                plt.subplot(234)
                plt.title('wklm1 (h2o)')
                plt.loglog(wklm1cols[:,col],pzmcols[1:,col],ls=linestyles[i1],label=str(fn))
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))
                if(legends_on==1):
                    plt.legend()

                plt.figure(i1+1)
                plt.subplot(235)
                plt.title('wklm2 (co2)')
                plt.semilogy(wklm2cols[:,col],pzmcols[1:,col],ls=linestyles[i1],label=str(fn))
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))
                if(legends_on==1):
                    plt.legend()

                plt.figure(i1+1)
                plt.subplot(236)
                plt.title('wklm3 (o3)')
                plt.semilogy(wklm3cols[1:,col],pzmcols[2:,col],ls=linestyles[i1],label=str(fn))
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))
                if(legends_on==1):
                    plt.legend()

                

                i3+=1

        i2 += 1

    # i1 += 1
    i_dir+=1


    # master indices: master[file,layer,column]
    tzm_master = np.array(tzm_master)
    pzm_master = np.array(pzm_master)
    altzm_master = np.array(altzm_master)
    boxlatcols_master = np.array(boxlatcols_master)
    lapsecritcols_master = np.array(lapsecritcols_master)
    meridtransp_master = np.array(meridtransp_master)
    abspncols_master = np.array(abspncols_master)
    A_oz_lcols_master=np.array(A_oz_lcols_master)
    abs_surf_lhcols_master=np.array(abs_surf_lhcols_master)
    totuflumcols_master=np.array(totuflumcols_master)
    abs_h2o_cols_master=np.array(abs_h2o_cols_master)
    abs_o3_cols_master=np.array(abs_o3_cols_master)
    abs_surf_cols_master=np.array(abs_surf_cols_master)
    d_mid_master=np.array(d_mid_master)
    d_trop_master=np.array(d_trop_master)
    rel_hum_master=np.array(rel_hum_master)
    wklm1_master=np.array(wklm1_master)
    wbrodlm_master=np.array(wbrodlm_master)


    box_abssw_tot_master = abs_surf_cols_master + abs_h2o_cols_master + abs_o3_cols_master
    boxtotnetflux_master = meridtransp_master[:,0,:] + box_abssw_tot_master[:,0,:] - totuflumcols_master[:,-1,:]

    # master indices for conv_trop_ind: master[file][column]
    conv_trop_ind_master = np.array(conv_trop_ind_master)


############################################################
print('Done')
# plt.tight_layout()
show()