    # Call the Fortran RCM Application from Python
# This version reads the MISR cloud data and other ERA-Interim data from .npy files created previously
# The .npy files have the maximum necessary dimensions, so the interpolation can happen on the fly (and quickly) in this script
# Thursday, 31 May 2018 

import os
import subprocess
import time
from subprocess import Popen, PIPE, STDOUT
import numpy as np
import math
from pylab import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import interpolate, stats
from scipy.interpolate import interp1d, interp2d, RectBivariateSpline, RegularGridInterpolator
from os import listdir
from time import localtime, strftime
from scipy import stats

ncols = 31
nlays = 30
days = 5000 #model days

def create_misr_cloud_inputs():
                
    cfs = misr_cf_latalt_max
    misr_alts = np.linspace(0,20,num=39) # altitudes in original MISR data
    misr_lats = np.linspace(-90,90,num=360)
    latbins = np.linspace(-90,90,num=ncols)
    altbin_edges = np.linspace(0,10,num=ncloudcols+1)
    altbins = np.zeros(ncloudcols)
    cld_taus = np.zeros(ncloudcols)

    for i in range(len(altbins)+1):
        altbins[i-1] = (altbin_edges[i-1] + altbin_edges[i])/2.0

    od_low = 3.0 * 0.0
    od_mid = 3.0 * 0.0
    od_high = 0.3 * 0.0

    for i in range(len(altbins)):
        if (altbins[i] < 3.0):
            cld_taus[i] = od_low
        elif(altbins[i] < 7.0):
            cld_taus[i] = od_mid
        else:
            cld_taus[i] = od_high

    binned_cf_lat = np.zeros((ncols,len(misr_alts)))
    cf_lat_bincounts = np.zeros((ncols,len(misr_alts)))
    binned_cf = np.zeros((ncols,ncloudcols))
    cf_bincounts = np.zeros((ncols,ncloudcols))

    for i in range(len(misr_lats)):
        ibin = np.argmin(abs(misr_lats[i] - collats))
        for j in range(len(misr_alts)):
            binned_cf_lat[ibin,j] = binned_cf_lat[ibin,j] + cfs[i,j]
            cf_lat_bincounts[ibin,j] = cf_lat_bincounts[ibin,j] + 1

    binned_cf_lat= binned_cf_lat / cf_lat_bincounts                
        
    for ibin in range(ncols):
        for j in range(len(misr_alts)):
            jbin = np.argmin(abs(misr_alts[j] - altbins))
            binned_cf[ibin,jbin] = binned_cf[ibin,jbin] + binned_cf_lat[ibin,j]
            cf_bincounts[ibin,jbin] = cf_bincounts[ibin,jbin] + 1

    clearfrac = np.zeros(ncols)

    binned_cf[0,0] = binned_cf[0,0] * 1.0

    for col in range(ncols):
        tempcloudfrac = 0.0
        for cloudcol in range(ncloudcols-1):
            tempcloudfrac = tempcloudfrac + binned_cf[col,cloudcol] - tempcloudfrac * binned_cf[col,cloudcol]
            clearfrac[col] = (1.0 - tempcloudfrac) / 2.0

    for col in range(ncols):
        filename = 'ccfracs col %2d' % (col)
        fileloc = outdir + filename
        file = open(fileloc,'w')
        for cloudcol in range(ncloudcols):
            file.write(str(binned_cf[col,cloudcol]))
            file.write('\n')
        file.write(str(clearfrac[col])) # add the clear fraction on at the end
        file.write('\n')
        file.close()

    for col in range(ncols):
        filename = 'ccalts col %2d' % (col)
        fileloc = outdir + filename
        file = open(fileloc,'w')
        for cloudcol in range(ncloudcols):
            file.write(str(altbins[cloudcol]))
            file.write('\n')
        file.write('1.0') # set altitude for clear column to 1 km
        file.write('\n')
        file.close()


    for col in range(ncols):
        filename = 'cctaus col %2d' % (col)
        fileloc = outdir + filename
        file = open(fileloc,'w')
        for cloudcol in range(ncloudcols):
            file.write(str(cld_taus[cloudcol]))
            file.write('\n')
        file.write('0.0') # set tau of clear column to 0
        file.write('\n')
        file.close()
def createlatdistbn(filename):
    fileloc = '/Users/nickedkins/Dropbox/Latitudinal Distributions/'+filename+'.txt'
    file = open(fileloc,'r')
    lat = []
    var = []
    with file as f:
        for l in f:
            lat.append(float(l.split(',')[0]))
            var.append(float(l.split(',')[1]))
    # lat = data[:,0]
    # var = data[:,1]    
    f = interpolate.interp1d(lat,var)    
    varinterp = list(f(latgrid))
    return varinterp
def findweightedglobavg(var):
    weightedglobalavgvar = sum(var*np.cos(radians(latgrid))) / sum(np.cos(radians(latgrid)))
    return weightedglobalavgvar
def interpolate_createprrtminput_lev(shortname,latparray,ps,lats):
    lats = lats
    pressures = ps
    xx,yy = np.meshgrid(lats[::-1],pressures)

    if(disttypelev[shortname] == 'lat'):

        z = latparray
        f = interpolate.RegularGridInterpolator((lats[::-1],pressures),z.T,bounds_error=False,fill_value=1000.0)
        xnew = latgrid
        ynew = pgrid
        xxnew, yynew = np.meshgrid(xnew,ynew)
        (xxnewr,yynewr) = (xxnew.ravel(),yynew.ravel())
        znew = f((xxnewr,yynewr),method='linear')
        znew=znew.reshape(nlays,ncols)
        znew = znew[:,::-1]

        if (shortname == pertvar):
            znew = znew * pert

            xnew = xnew[::-1]
            ynew = ynew[::-1]

    elif(disttypelev[shortname]=='avg'):
            
        z = latparray
        f = interpolate.RegularGridInterpolator((lats[::-1],pressures),z.T,bounds_error=False,fill_value=1000.0)

        xnew = latgrid
        ynew = pgrid
        xxnew, yynew = np.meshgrid(xnew,ynew)
        (xxnewr,yynewr) = (xxnew.ravel(),yynew.ravel())
        znew = f((xxnewr,yynewr),method='linear')
        znew=znew.reshape(nlays,ncols)

        znew = znew[:,::-1]

        if (shortname == pertvar):
            znew = znew * pert


            xnew = xnew[::-1]
            ynew = ynew[::-1]

            zavg = np.zeros(nlays)

        for col in range(0,ncols):
            zavg = zavg + znew[:,col] * latweights[col]    

        zavg = zavg / sum(latweights)


    #f = interpolate.interp1d(pressures,zavg)        
    #znew = f(pgrid)

    if (disttypelev[shortname]=='lat'):            

        # Write input files for PRRTM
        for col in range(ncols):

            filename = '%s vert col %2d' % (shortname,col)
            fileloc = outdir + filename
            file = open(fileloc,'w')

            for i in range(len(znew[:,col])):
                file.write(str(znew[i,col]))
                file.write('\n')

            file.write('&')

            file.close()

    elif (disttypelev[shortname]=='avg'):

        for col in range(ncols):

            filename = '%s vert col %2d' % (shortname,col)
            fileloc = outdir + filename
            file = open(fileloc,'w')

            for i in range(nlays):
                file.write(str(zavg[i]))
                file.write('\n')

            file.close()
def interpolate_createprrtminput_sfc(shortname,latarray,lats):
    lats = lats
    z = latarray
    f = interp1d(lats,z)
    varss_int = f(latgrid)

    if (shortname == pertvar):
        varss_int = varss_int * pert

    if (disttypesfc[shortname]=='avg'):
        zavg = 0.0
        for col in range(0,ncols):
            zavg = zavg + varss_int[col] * latweights[col]        
            zavg = zavg / sum(latweights)
            filename = '%s lats' % (shortname)
            fileloc = outdir + filename
            file = open(fileloc,'w')
            for i in range(len(varss_int)):
                file.write(str(zavg))
                file.write('\n')
            file.close()

    elif(disttypesfc[shortname]=='lat'):
        filename = '%s lats' % (shortname)
        fileloc = outdir + filename
        file = open(fileloc,'w')

        for i in range(len(varss_int)):
            file.write(str(varss_int[i]))
            file.write('\n')

        file.close()
def create_manual_cloud_inputs():
    
    # Calculate clear sky fraction with random overlap assumption
    c1 = 0.
    c2 = 0.
    for cloudcol in range(ncloudcols):
        c1 = manual_clouds[cloudcol][1]
        ctot = c1 + c2 - c1*c2
        c2 = ctot
    clearfrac = np.ones(ncols) * ( 1.0 - ctot )
    manual_clouds_wghtd_frac = np.zeros(ncloudcols)

    for cloudcol in range(ncloudcols):
        manual_clouds_wghtd_frac[cloudcol] = manual_clouds[cloudcol][1] * ctot / sum(np.array(manual_clouds)[:,1])

    for col in range(ncols):
        filename = 'ccfracs col %2d' % (col)
        fileloc = outdir + filename
        file = open(fileloc,'w')
        for cloudcol in range(ncloudcols):
            file.write( str( manual_clouds_wghtd_frac[cloudcol] ) )
            file.write('\n')
        file.write(str(clearfrac[col])) # add the clear fraction on at the end
        file.write('\n')
        file.close()

    for col in range(ncols):
        filename = 'ccalts col %2d' % (col)
        fileloc = outdir + filename
        file = open(fileloc,'w')
        for cloudcol in range(ncloudcols):
            file.write( str( manual_clouds[cloudcol][0] ) )
            file.write('\n')
        file.write('1.0') # set altitude for clear column to 1 km
        file.write('\n')
        file.close()


    for col in range(ncols):
        filename = 'cctaus col %2d' % (col)
        fileloc = outdir + filename
        file = open(fileloc,'w')
        for cloudcol in range(ncloudcols):
            file.write( str( manual_clouds[cloudcol][2] ) )
            file.write('\n')
        file.write('0.0') # set tau of clear column to 0
        file.write('\n')
        file.close()

#project_dir = '/Users/nickedkins/Dropbox/GitHub Repositories/2D-RCM-Home/2D-RCM/'
project_dir = '/Users/nickedkins/Dropbox/GitHub Repositories/2D-RCM/2D-RCM/'

interpdir = '/Users/nickedkins/Dropbox/Input Data for RCM Interpolation/'
outdir = project_dir+'Input Distributions/' #output file directory

misr_cf_latalt_max = np.load(interpdir+'misr_cf_latalt.npy') # maximum-sized array of cloud data, to interpolate onto new smaller grid here
misr_cf_latalt_max = np.nan_to_num(misr_cf_latalt_max)
xbounds = np.linspace(-1,1,ncols+1)
latbounds = np.rad2deg(np.arcsin(xbounds))
collats = np.zeros(ncols)
colwidth = 180.0 / (ncols+1)

for i in range(1,len(latbounds)):
    collats[i-1] = (latbounds[i-1] + latbounds[i])/2

latweights = np.cos(radians(collats))
latgrid = collats
pgrid = np.linspace(1000,1,nlays)

q_latp_max = np.load(interpdir+'q_latp.npy')
o3_latp_max = np.load(interpdir+'o3_latp.npy')
fal_lat_max = np.load(interpdir+'fal_lat.npy')

q_ps = np.load(interpdir+'q_ps.npy')
o3_ps = np.load(interpdir+'o3_ps.npy')

q_lats = np.load(interpdir+'q_lats.npy')
o3_lats = np.load(interpdir+'o3_lats.npy')
fal_lats = np.load(interpdir+'fal_lats.npy')

pa = 0.3    
sc = [1362.0]

#pico2s = np.linspace(400e-6,3200e-6,num=5)

#pico2s = np.logspace(-4,2,num=5,base=10.0)
pico2s = [400e-6]
# pin2s = np.logspace(-1,2,num=10,base=10.0)
pin2s = [1.0]

#pico2s = [400e-6,3200e-6]


pertlay=0
pertcol=0
pertvar = 'fal'
pert =1.0

#inversion_strengths = np.linspace(0,-40,num=10)
inversion_strength = [0.0]
inversion_col = 1
#polar_tgs = [240,250,260,270,280]

# sebfacs = np.linspace(0.0,1e-1,num=5)
#sebfacs = [0.2]


surf_layer_depth = 1e1
surf_layer_density = 3.e3
surf_layer_shc = 3.e3

sebfac = (60.*60.*24.) / (surf_layer_depth * surf_layer_density * surf_layer_shc)

# sas = np.linspace(0.2,0.8,num=5)
sas = [0.4]
#tboundms = np.linspace(250,335,num=10)
tboundms = [288.0]

#for pertcol in range(ncols):

# global_lapses = np.linspace(-2,-8,5)
global_lapses = [-5.7]

cloud_loc_type = 0 # 0: pressure (hPa), 1: altitude (km), 2: temperature (K)


cloud_source = 0

cld_heights = np.linspace(1,10,10)
cld_heights = [5.0]

for cld_height in cld_heights:

    manual_clouds = []
    manual_clouds.append([cld_height,0.5,3.0])
    ncloudcols = shape(manual_clouds)[0]

    for tboundm in tboundms:
    
        for sa in sas:
    
            sa = [sa] * ncols
            
    
            for pin2 in pin2s:
    
                for pico2 in pico2s:

                    if ( cloud_source == 0 ):

                        create_manual_cloud_inputs()

                    elif ( cloud_source == 1 ):

                        create_misr_cloud_inputs()
                
                    #pertvars = ['q','cc','ciwc','clwc','o3']
                    pertvars = ['q']
    
                    # shortnameslev = ['cc','clwc','o3','q','ciwc']
                    shortnameslev = ['q', 'o3']
                    longnameslev = {'cc':'Cloud fraction','clwc':'Cloud liquid water content (kg/kg)','o3':'Ozone mixing ratio','q':'Specific humidity (kg/kg)','ciwc':'Cloud ice water content (kg/kg)'}
                    #disttypelev = {'cc':'lat','clwc':'lat','o3':'lat','q':'lat','ciwc':'lat'}
                    #disttypelev = {'cc':'lat','clwc':'lat','o3':'lat','q':'lat','ciwc':'lat'}
                    disttypelev = {'cc':'lat','clwc':'lat','o3':'lat','q':'lat','ciwc':'lat'}
                
                    shortnamessfc = ['fal']
                    longnamessfc = {'fal':'Surface albedo'}
                    disttypesfc = {'fal':'lat'}
                
                    loop = 1
                
                    print("Creating input files by interpolating ERA-Interim data")
                
                    interpolate_createprrtminput_lev('q',q_latp_max,q_ps,q_lats)
                    interpolate_createprrtminput_lev('o3',o3_latp_max,o3_ps,o3_lats)
                    interpolate_createprrtminput_sfc('fal',fal_lat_max,fal_lats)
                
                    lc = createlatdistbn('Doug Mason Lapse Rate vs Latitude')
                    #lc[0] = -10.0
                    #lc[ncols-1] = -10.0
                    #for i in range(len(lc)):
                    #    lc[i] = lcmean

                    lch = createlatdistbn('Cloud Top Height')
                    srh = createlatdistbn('Relative Humidity')
                    # srh = [0.8] * ncols
                    sa = createlatdistbn('Surface Reflectance')
                    #sa = [0.5] * ncols
                    lcf = createlatdistbn('Cloud Fraction')
                    lcod = createlatdistbn('Cloud Optical Thickness')
                    tg = createlatdistbn('Surface Temperature')
                    # tg = [tboundm] * ncols
                    # tg = [250.] * ncols
                
                    #lc = [-5.88]*ncols
                    #lch = [4.46]*ncols
                    #srh = [0.787]*ncols
                    #lcf = [0.657]*ncols
                    #lcod = [3.978]*ncols
                
                    mc = pico2  
                
                    ur = 0.5
                    icldm = 1
                    rmin = 3e-6
                    hct = 230.0
                    hcf = 0.04e-9   
                    hcod = 0.7e-9
                    mct = 235.0
                    mcf = 0.05e-9
                    mcod = 1.5e-9
                    #lct = 250.0
                    #lcf = 0.5
                    #lcod = 5.0
                    tp = 5.0
                    #fth = np.zeros(ncols)
                    #for i in range(ncols):
                    #    fth[i] = 15.0 - abs(collats[i])/18.0
                    fth = [500.] * ncols
                    ol = nlays
                    asp = 2.0   
                    cs = 0
                    pbo = 0 
                    fswon = 0  
                    fsw = 239.4
                    fp = 0
                    ps1 = 0
                    af = 1.0
                    dalr = 0 #convection type
                    npb = 1
                    o3sw = 1
                    h2osw = 0
                    nl = nlays
                    maxhtr = 0.1
                    asf = 4.0
                    tuf = 1.0   
                    n2inv = pin2
                    #n2inv = 0.8
                    o2inv = 0.0
                    htransp = 1.0 #reduce lapse rate to account for horizontal transport
                    ipe = 1
                    dp = 1
                    mtranspfac = 2.0 * 0.0
                    boxnetfluxfac = 0.2
                    twarm = 288
                    tcold = 268
                    phim = 45 * 3.14 / 180
                    ks = 0.25
                    kl = 0.25
                    eta = 0.75
                    planet_radius = 6.37e6
                    planet_rotation = 7.29e-5
                    t_min = 100.0
                    sfc_heating = 0 #surface energy budget warms/cools surface? 1=yes, 0=no
                    playtype = 0 #pressure layer type. 0=equal p thickness, 1=sigma
                    ur_htr = 0.5
                    ur_toafnet = 50.
                    ur_seb = 1e10
                    couple_tgta = 1
                
                    ur1 = ur
                
                    counter = 0
                
                    ur = ur1
                
                    params = [ncols,ncloudcols+1,pa,sc,tg,lc,days,mc,ur,icldm,rmin,hct,hcf,hcod,mct,mcf,mcod,lch,lcf,lcod,tp,sa,list(fth),ol,asp,cs,pbo,fswon,fsw,fp,srh,ps1,af,dalr,
                    npb,o3sw,h2osw, nl, maxhtr, asf, tuf, pico2, n2inv, o2inv, htransp, ipe, dp, mtranspfac,boxnetfluxfac,pertlay,pertcol,list(collats),inversion_strength,inversion_col,
                    twarm,tcold,phim,ks,kl,eta,planet_radius,planet_rotation,list(latbounds),t_min,sebfac,sfc_heating,playtype,ur_htr,ur_toafnet,ur_seb,couple_tgta]
                    
                    f = open(project_dir+'/Earth RCM Parameters','w')
                    for m in params:
                        if type(m) is list:
                            for i in range(len(m)-1):
                                f.write(str(m[i]))
                                f.write(',')
                            f.write(str(m[len(m)-1]))
                            f.write('\n')
                        else:
                            f.write(str(m))
                            f.write('\n')
                
                    f.write('$')
                    f.close()
                
                    print("Calling PRRTM")
                
                    time.sleep(2)
                
                    for i in range(1,2):
                
                        loc = project_dir+'2D RCM GitHub'
                        os.chdir(project_dir)
                        print(os.getcwd())  # Prints the current working directory
                        p = subprocess.Popen([loc])
                
                        stdoutdata, stderrdata = p.communicate()
                        # print 'return code = %4d' % (p.returncode)
                        print('return code = {}'.format(p.returncode))
                        print('------------------------------------------------------------------------------------------')
                        print
                
                        if (p.returncode == 0):
                            break
                        elif (p.returncode == -11 or p.returncode == -8 or p.returncode == 11 or p.returncode == 12):
                            ur = ur*2
                        continue
                
                        counter = counter + 1

########################################################################################################################

os.system('say "Done"')
show()