# Read and plot 2D RCM output data

# Just plot the ouput for the current run of PRRTM
import numpy as np
import matplotlib as matplotlib
import matplotlib.pyplot as plt
from pylab import *
from scipy import interpolate
from os import listdir
# import pandas as pd

plot_all_vert_profiles = 0
legends_on = 0
grids_on = 1

directories = [
'_Current Output/'
]

directories = [
'/Users/nickedkins/Dropbox/GitHub Repositories/Uni/2D-RCM/_Useful Data/mtransp expts/h2o=erai/latfac/'
]

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
    print "x:", x
    print "lats:", lats
    print "sinlats:", np.cos(np.deg2rad(lats))
    print "sumlats:", sum(np.cos(np.deg2rad(lats)))
    x_avg = sum(x * np.cos(np.deg2rad(lats))) / sum(np.cos(np.deg2rad(lats)))
    return x_avg


obs_file = '/Users/nickedkins/Dropbox/Figure Data/shinesinha_deltaTvsp.txt'
obs_data = np.genfromtxt(obs_file,delimiter=',')

t_obs = obs_data[:,0]
p_obs = obs_data[:,1]

p_obs, t_obs = zip(*sorted(zip(p_obs, t_obs)))

# plt.figure(1)
# plt.plot(t_obs,p_obs,'--',label='Shine and Sinha')
# plt.ylim(1000,0)



# plt.figure(1)
# plt.subplot(341)
# plt.plot(t_obs,z_obs,'--',label='ERA-Interim')
# plt.plot(t_grey,z_grey,label='Grey')
# plt.xlabel('Temperature (K)')
# plt.ylabel('Altitude (km)')
# plt.ylim(0,20)

# plt.figure(1)
# plt.subplot(222)
# plt.title('lapse')
# dt = np.zeros(len(t_obs))
# dz = np.zeros(len(z_obs))



# for i in range(3,len(t_obs)):
#     dt[i] = t_obs[i] - t_obs[i-3]
#     dz[i] = z_obs[i] - z_obs[i-3]
# plt.semilogy(dt/dz,p_obs)
# plt.ylim(max(p_obs),min(p_obs))

# plt.semilogy(t_obs,p_obs,'--')
# plt.ylim(max(p_obs),min(p_obs))


linestyles = ['-','--','--','o']

colors = ['b','r','g','orange','purple','yellow','pink']



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

    sol_inc = 1362.0/4.0
    # abs_h2o = sum(abspncols[:nlayersm-1,:])*sol_inc/ncols / factor
    # abs_o3 = sum(A_oz_lcols[2:nlayersm,:])*sol_inc/ncols / factor
    # abs_surf = np.mean(abs_surf_lhcols[0,:]) / factor


    return tzmcols,pzmcols,wklm1cols,totuflumcols,htrmcols,altzmcols,pavelmcols,htro3cols,totdflumcols,wklm2cols,A_oz_lcols,\
    abspncols,abs_surf_lhcols,tboundmcols,tavelmcols,nlayersm,ncols,boxlatcols,htrh2ocols,wklm3cols,convcols,wbrodlmcols,\
    lapsecritcols,meridtransp,abs_h2o_cols,abs_o3_cols,abs_surf_cols

# nlayersms=[31,100]
# ncols=5


i1 = 0

tzm_store = np.zeros( ( 200, len(directories), 100 ) )
pico2_store = np.zeros( (len(directories), 100 ) )

cld_heights = np.linspace(1,10,10)

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

    filenames = []

    dir_label = directory.split('/')[-2]
    print dir_label

    ls = linestyles[i1]
    color = colors[i1]

    a = sorted(listdir(directory))

    filenames.append(a)

    counter=0

    i2 = 0

    for fn in a:
        if (fn=='.DS_Store'):
            continue
        tzmcols,pzmcols,wklm1cols,totuflumcols,htrmcols,altzmcols,pavelmcols,htro3cols,totdflumcols,wklm2cols,A_oz_lcols,abspncols,\
        abs_surf_lhcols,tboundmcols,tavelmcols,nlayersm,ncols,boxlatcols,htrh2ocols,wklm3cols,convcols,wbrodlmcols,lapsecritcols,\
        meridtransp,abs_h2o_cols,abs_o3_cols,abs_surf_cols = readfile(fn,counter)

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

                print p_trop,',', t_trop,',', z_trop/1000.,',', pzmcols[0,col],',', tzmcols[0,col]

                # plt.figure(2)
                # plt.plot(lapsecritcols[0,col],z_trop/1000.,'o')

                # for i in range(len(altzmcols[:,col])):
                #     print altzmcols[i,col]/1000., ',', pzmcols[i,col], ',', tzmcols[i,col]

                plt.figure(i1+1)
                # plt.figure(1)
                plt.subplot(341)
                plt.title('tzm')
                # plt.plot(tzmcols[:,col],altzmcols[:,col]/1000.,ls=linestyles[i1],label='Spectral '+str(fn))
                plt.semilogy(tzmcols[:,col],pzmcols[:,col],label=str(fn))
                # plt.plot(tzmcols[:,col],altzmcols[:,col],'-o',label=str(fn))
                # plt.plot(tzmcols[conv_trop_ind,col],altzmcols[conv_trop_ind,col]/1000.,'*',markersize=20)
                plt.plot(tzmcols[conv_trop_ind,col],pzmcols[conv_trop_ind,col],'*',markersize=20)
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))
                if(grids_on==2):
                    plt.gca().minorticks_on()
                if(grids_on==1):
                    plt.grid(which='both',axis='both')
                if(legends_on==1):
                    plt.legend()    

                # dt_obs = np.zeros(nlayersm)
                # dz_obs = np.zeros(nlayersm)
                # plt.subplot(222)
                # for i in range(1,nlayersm):
                #     dt_obs[i] = tzmcols[i,col] - tzmcols[i-1,col]
                #     dz_obs[i] = (altzmcols[i,col] - altzmcols[i-1,col]) / 1000.
                # plt.plot(-dt_obs/dz_obs,pzmcols[1:,col])


                
                plt.figure(i1+1)
                plt.subplot(342)
                plt.title('htrm')
                plt.semilogy(htrmcols[:,col],pzmcols[:,col],ls=linestyles[i1],label=str(fn))
                plt.xlim(-5,5)
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))
                plt.axvline(-0.1,ls='--')
                plt.axvline(0.1,ls='--')
                if(legends_on==1):
                    plt.legend()
                
                plt.figure(i1+1)
                plt.subplot(343)
                plt.title('totuflum')
                plt.semilogy(totuflumcols[:,col],pzmcols[:,col],ls=linestyles[i1],label='up '+dir_label)
                # plt.semilogy(totdflumcols[:,col],pzmcols[:,col],ls=linestyles[i1],label='down '+dir_label)
                # plt.semilogy(totdflumcols[:,col]-totuflumcols[:,col],pzmcols[:,col],ls=linestyles[i1],label='net '+dir_label)
                plt.axvline(0,ls='--')
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))
                if(legends_on==1):
                    plt.legend()

                plt.figure(i1+1)
                plt.subplot(345)
                plt.title('abs_o3')
                plt.semilogy(A_oz_lcols[:,col]*1362./4.,pzmcols[1:,col],ls=linestyles[i1])
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))
                if(legends_on==1):
                    plt.legend()

                plt.figure(i1+1)
                plt.subplot(3,4,7)
                plt.title('wklm1 (h2o)')
                plt.loglog(wklm1cols[:,col],pzmcols[1:,col],ls=linestyles[i1],label=str(fn))
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))
                if(legends_on==1):
                    plt.legend()

                plt.figure(i1+1)
                plt.subplot(3,4,8)
                plt.title('wklm2 (co2)')
                plt.semilogy(wklm2cols[:,col],pzmcols[1:,col],ls=linestyles[i1],label=str(fn))
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))
                if(legends_on==1):
                    plt.legend()

                plt.figure(i1+1)
                plt.subplot(3,4,9)
                plt.title('wklm3 (o3)')
                plt.semilogy(wklm3cols[1:,col],pzmcols[2:,col],ls=linestyles[i1],label=str(fn))
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))
                if(legends_on==1):
                    plt.legend()

                plt.figure(i1+1)
                plt.subplot(3,4,10)
                plt.title('wbrodlm')
                plt.semilogy(wbrodlmcols[:,col],pzmcols[1:,col],ls=linestyles[i1])
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))
                if(legends_on==1):
                    plt.legend()

                plt.figure(i1+1)
                plt.subplot(3,4,11)
                plt.title('htro3')
                plt.semilogy(htro3cols[:,col],pzmcols[1:,col],ls=linestyles[i1],label=str(fn))
                plt.xlim(-5,5)
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))
                if(legends_on==1):
                    plt.legend()

                plt.figure(i1+1)
                plt.subplot(3,4,12)
                plt.title('htrh2o')
                plt.semilogy(htrh2ocols[:,col],pzmcols[1:,col],ls=linestyles[i1],label=str(fn))
                # plt.xlim(-5,5)
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))
                if(legends_on==1):
                    plt.legend()

                

                plt.figure(i1+1)
                plt.subplot(3,4,6)
                plt.title('htrmlw')
                plt.semilogy(htrmlwcols[:,col],pzmcols[1:,col],ls=linestyles[i1])
                plt.xlim(-5,5)
                plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))
                plt.axvline(-0.03,ls='--')
                plt.axvline(0.03,ls='--')
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


    # master indices for conv_trop_ind: master[file][column]
    conv_trop_ind_master = np.array(conv_trop_ind_master)



    # i_fn = 0
    # for fn in a:
    #     if (fn=='.DS_Store'):
    #             continue

    #     plt.figure(1)
    #     plt.subplot(221)
    #     plt.plot(boxlatcols_master[i_fn,0,:],-1.0*lapsecritcols_master[i_fn,0,:],label=str(fn))
    #     plt.xlabel('Latitude')
    #     plt.ylabel('Lapse rate (K/km)')
    #     plt.legend()

    #     plt.subplot(222)
    #     plt.plot(boxlatcols_master[i_fn,0,:],tzm_master[i_fn,conv_trop_ind_master[i_fn,:],range(ncols)],'-',label=str(fn))
    #     plt.xlabel('Latitude')
    #     plt.ylabel('Tropopause temperature (K)')
    #     plt.legend()
        
    #     plt.subplot(223)
    #     plt.plot(boxlatcols_master[i_fn,0,:],altzm_master[i_fn,conv_trop_ind_master[i_fn,:],range(ncols)]/1000.,'-',label=str(fn))
    #     plt.xlabel('Latitude')
    #     plt.ylabel('Tropopause altitude (km)')
    #     plt.legend()

    #     # plt.subplot(224)
    #     # plt.plot(boxlatcols_master[i_fn,0,:],pzm_master[i_fn,conv_trop_ind_master[i_fn,:],range(ncols)],'-',label=str(fn))
    #     # plt.xlabel('Latitude')
    #     # plt.ylabel('Tropopause pressure (hPa)')
    #     # plt.legend()

    #     plt.subplot(224)
    #     plt.plot(boxlatcols_master[i_fn,0,:],tzm_master[i_fn,0,range(ncols)],'-',label=str(fn))
    #     plt.xlabel('Latitude')
    #     plt.ylabel('Surface temperature (K)')
    #     plt.legend()    


    #     i_fn+=1




    # filenames = np.array(filenames[0])

    # x = np.linspace(-90,90,ncols)

    lats = boxlatcols_master[:,0,:]
    tgs = tzm_master[:,0,:]

    abs_sw_tot = abs_h2o_cols_master + abs_o3_cols_master + abs_surf_cols_master

    # print latwghtavg()

    # for i_fn in range(len(a)-1):
    #     plt.plot(boxlatcols_master[i_fn,0,:],tzm_master[i_fn,0,:])
    #     print latwghtavg(tzm_master[i_fn,0,:],boxlatcols_master[i_fn,0,:])

    for i_fn in range(0,len(boxlatcols_master[:,0,0])):

        plt.figure(1)
        plt.subplot(121)
        plt.plot(boxlatcols_master[i_fn,0,:],meridtransp_master[i_fn,0,:],'-o',label=a[i_fn])
        plt.xlabel('Latitude')
        plt.ylabel('$\Delta$ box meridional transport (Wm$^{-2}$)')
        plt.axhline(0)
        plt.legend()

        plt.subplot(122)
        plt.plot(boxlatcols_master[i_fn,0,:],tzm_master[i_fn,0,:],'-o',label=a[i_fn])
        plt.xlabel('Latitude')
        plt.ylabel('T_{g}')
        plt.legend()


    



    # plt.subplot(221)
    # plt.plot(boxlatcols_master[0,0,:],tzm_master[1,0,:] - tzm_master[0,0,:],'-o')
    # plt.axhline(0)
    # plt.xlabel('Latitude')
    # plt.ylabel('$\Delta T_g$')

    # plt.subplot(222)
    # plt.plot(boxlatcols_master[0,0,:],meridtransp_master[1,0,:] - meridtransp_master[0,0,:],'-o')
    # plt.axhline(0)
    # plt.xlabel('Latitude')
    # plt.ylabel('$\Delta mtransp$')



    # print tzm_master[:,0,:], boxlatcols_master[0,0,:], latwghtavg(tzm_master[:,0,:],boxlatcols_master[0,0,:])

    

    # plt.figure(1)
    # plt.subplot(221)
    # plt.plot(lats,tzm_master[0,0,range(ncols)])
    # plt.plot(lats,tzm_master[1,0,range(ncols)])
    # plt.subplot(222)
    # plt.plot(lats,tzm_master[0,conv_trop_ind_master[0,:],range(ncols)])
    # plt.plot(lats,tzm_master[1,conv_trop_ind_master[1,:],range(ncols)])
    # plt.legend()

    # pperts = np.linspace(1000,0,10)

    # if (dir_label == 'Manabe-Wetherald 3 Clouds' or dir_label == 'Manabe-Wetherald Warmer'):
    #     pperts = np.linspace(1000,0,10)

    # # plt.figure(1)
    # # for i_fn in range(len(a)-1):
    # #     if(a[i_fn]=='.DS_Store'):
    # #         continue
    # plt.plot(tzm_master[1:,0,:]-tzm_master[0,0,:],pperts,label=dir_label)
    # plt.xlabel('$\Delta T_{surf}$ (K)')
    # plt.ylabel('Pressure at bottom of 50 hPa H$_2$O perturbation region (hPa)')
    # plt.ylim(1000,0)

    # print tzm_master[:,0,:]


plt.legend()

# for i in range( shape(tzm_master)[0] ):

#     lats = boxlatcols_master[i,0,:]
#     pzms = pzm_master[i,:,0]
#     tzms = tzm_master[i,:,:]
#     altzms = altzm_master[i,:,0]
#     # tsgm_weighted = sum(tzms[0,:] * np.cos(np.deg2rad(lats)) / sum(np.cos(np.deg2rad(lats)) ))



    # plt.figure(1)
    # plt.subplot(121+i)
    # plt.contourf(lats,pzms,tzms,20)
    # plt.gca().set_yscale('log')
    # plt.ylim(1000,10)
    # plt.xlabel('Latitude')
    # plt.ylabel('Altitude')
    # plt.colorbar()

    # plt.figure(2)
    # plt.plot(lats,tzms[0,:])
    

# plt.subplot(223)
# plt.gca().set_yscale('log')
# plt.ylim(1000,10)
# plt.contourf(x,y,tzm_master[1,:,:]-tzm_master[0,:,:],20)
# plt.xlabel('Latitude')
# plt.ylabel('Altitude')
# plt.colorbar()



############################################################
print 'Done'
plt.tight_layout()
show()