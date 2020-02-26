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

# project_dir = '/Users/nickedkins/Dropbox/GitHub Repositories/Uni/2D-RCM/'
project_dir = '/Users/nickedkins/Dropbox/GitHub Repositories/Home/2D-RCM/'

os.chdir(project_dir)

ncolss = [7]
ncloudcolss = [2]
nlayss = [600]
od_low = 3.0
od_mid = 3.0
od_high = 3.0
nperts = 1
timesteps = 3000
ur_htr = 0.3
# ur_htr = 3.0
days = timesteps/ur_htr
min_press = 1.
cloud_source = 0 #0 for manual, 1 for MISR
steps_before_first_eqbcheck = 200
snapshot=0
h2o_sources=[1] # 0=ERA-I mixh2o, 1=MW67 RH, 2=Cess RH, 3=Kasting, 4=Ramirez, 5=constant with lat, 6=Kluft19
o3_sources = [1]	 #1=erai, 2=RCEMIP
maxhtr = 0.03
lcs = [-5.7]
lapse_types = [2] # 0=critical lapse rate, 1=H82, 2=Mason
convecttype = 0 #convection type. 0: normal critical lapse 1: for forcing expt, convect to fixed ptrop and no higher 2: MALR
manual_clouds = []
# manual_clouds.append([750.,0.5,3.])
fswon = 0 #0:off, 1:on
fsws = [240.] #238.24 to replicate RD
mtransp_types = [2] #1=simple diffusion, 2=Vladilo
# gas_amt_fac_co2s = [1/2.,1.,2.,4.,8.,16.]
gas_amt_fac_co2s = [1.]
gas_amt_fac_ch4s = [1.]		
gas_amt_fac_h2os = [1.]
forcing_expt = 0 #0=normal, 1=forcing expt: stratosphere adjusts, surface and tropopshere are fixed, ptrop fixed
tp = 0.1
if (forcing_expt==1):
	tp = tp * 1e12



for nlays in nlayss:
	nlays = int(nlays)
	# ur_htr = nlays/1000.
	ur_htr = nlays/500.

	for ncols in ncolss:
		ncols = int(ncols)

		for ncloudcols in ncloudcolss:

			def create_misr_cloud_inputs(add_cld_alt):
							
				cfs = misr_cf_latalt_max
				misr_alts = np.linspace(0,20,num=39) # altitudes in original MISR data
				misr_lats = np.linspace(-90,90,num=360)
				latbins = np.linspace(-90,90,num=ncols)
				altbin_edges = np.linspace(0,20,num=ncloudcols+1)
				altbins = np.zeros(ncloudcols)
				cld_taus = np.zeros(ncloudcols)

				for i in range(len(altbins)+1):
					altbins[i-1] = (altbin_edges[i-1] + altbin_edges[i])/2.0



				for i in range(len(altbins)):
					if (altbins[i] < 3.0):
						cld_taus[i] = od_low
					elif(altbins[i] < 10.):
						cld_taus[i] = od_mid
					else:
						cld_taus[i] = od_high

				print(altbins, cld_taus)

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
					for cloudcol in range(ncloudcols):
						tempcloudfrac = tempcloudfrac + binned_cf[col,cloudcol] - tempcloudfrac * binned_cf[col,cloudcol]
						# clearfrac[col] = (1.0 - tempcloudfrac) / 2.0
						# clearfrac[col] = (1.0 - tempcloudfrac)

					clearfrac[col] = 1.0 - np.sum(binned_cf[col,:])

				for col in range(ncols):
					filename = 'ccfracs col %2d' % (col)
					fileloc = outdir + filename
					file = open(fileloc,'w')
					for cloudcol in range(ncloudcols):
						file.write(str(binned_cf[col,cloudcol]))
						file.write('\n')
					file.write('0.1')
					file.write('\n')
					file.write(str(clearfrac[col])) # add the clear fraction on at the end
					file.write('\n')
					file.close()

				for col in range(ncols):
					filename = 'ccalts col %2d' % (col)
					fileloc = outdir + filename
					file = open(fileloc,'w')
					for cloudcol in range(ncloudcols):
						file.write(str(altbins[cloudcol]+add_cld_alt)) #nje
						file.write('\n')
					file.write('15.0')
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
					file.write('0.3')
					file.write('\n')
					file.write('0.0') # set tau of clear column to 0
					file.write('\n')
					file.close()
			def createlatdistbn(filename):
				fileloc = 'Latitudinal Distributions/'+filename+'.txt'
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
				templats = np.linspace(-90,90,30)
				varinterp = np.zeros(len(latgrid))
				cossums=np.zeros(len(latgrid))
				for j in range(len(latgridbounds)-1):
					for i in range(len(templats)):
						if (templats[i] > latgridbounds[j] and templats[i] < latgridbounds[j+1]):
							varinterp[j] += f(templats[i]) * np.cos(np.deg2rad(templats[i]) )
							cossums[j] += np.cos(np.deg2rad(templats[i]))
				varinterp = list(varinterp / cossums)
				return varinterp
			def findweightedglobavg(var):
				weightedglobalavgvar = sum(var*np.cos(radians(latgrid))) / sum(np.cos(radians(latgrid)))
				return weightedglobalavgvar
			def interpolate_createprrtminput_lev(shortname,latparray,ps,lats,lat_facs):
				lats = lats #latitudes at 5 deg spacing - resolution at which I sampled ERA-I earlier 
				pressures = ps # ERA-I pressures
				xx,yy = np.meshgrid(lats[::-1],pressures) # grid of lat-p in my sampled ERA-I resn

				if(disttypelev[shortname] == 'lat'):

					z = latparray # values of the variable on the maximum resn lat-p grid that I sampled from ERA-I earlier
					f = interpolate.RegularGridInterpolator((lats[::-1],pressures),z.T,bounds_error=False,fill_value=1000.0) # function that will return variable on whatever latp grid is passed to it
					# This is where I need to integrate instead of evaluating at a point
					xnew = latgridbounds # latgrid is the lats at the centre of each lat col in the 2D RCM. Need the edges
					ynew = pgrid # pgrid is the ps in the 2D RCM
					# Think/read about how to bin the large array into the small one.
					znew = np.zeros( (len(latgridbounds)-1, len(pgrid)-1) )
					weights = np.zeros( (len(latgridbounds)-1, len(pgrid)-1) )
					sinlat_int = np.linspace(-1.,1.,ncols*2)
					# lats_int = np.rad2deg(np.arcsin(sinlat_int))
					lats_int = np.linspace(-90,90,ncols*2) # lats to integrate over (step size)
					# lats_int = [0] # lats to integrate over (step size)
					pressures_int = np.linspace(1000,1,nlays*2) # ps to integrate over (step size)
					for i_lat in range(len(lats_int)):
						for i_p in range(len(pressures_int)):
							for i_latg in range(len(latgridbounds)-1):
								for i_pg in range(len(pgrid)-1):
									if (latgridbounds[i_latg] <= lats_int[i_lat] < latgridbounds[i_latg+1]):
										if (pgrid[i_pg] >= pressures_int[i_p] > pgrid[i_pg+1]):
											znew[i_latg,i_pg] += f((lats_int[i_lat],pressures_int[i_p]),method="linear") * np.cos(np.deg2rad(lats_int[i_lat]))
											weights[i_latg,i_pg] += np.cos(np.deg2rad(lats_int[i_lat]))

					znew = znew/weights


					# xxnew, yynew = np.meshgrid(xnew,ynew)
					# (xxnewr,yynewr) = (xxnew.ravel(),yynew.ravel())
					# znew = f((xxnewr,yynewr),method='linear')
					# znew=znew.reshape(nlays,ncols)
					# znew = znew[:,::-1]

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

						for i in range(len(znew[col,:])):
							file.write(str(znew[col,i]*lat_facs[col]))
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
			def interpolate_createprrtminput_sfc(shortname,latarray,lats,lat_facs):
				lats = lats
				z = latarray
				f = interp1d(lats,z)

				weights = np.zeros( len(latgridbounds)-1 )
				lats_int = np.linspace(-90,90,ncols*2.)
				varss_int = np.zeros(len(latgridbounds)-1)
				for i_lat in range(len(lats_int)):
					for i_latg in range(len(latgridbounds)-1):
						if (latgridbounds[i_latg] <= lats_int[i_lat] < latgridbounds[i_latg+1]):
							varss_int[i_latg] += f(lats_int[i_lat]) * np.cos(np.deg2rad(lats_int[i_lat]))
							weights[i_latg] += np.cos(np.deg2rad(lats_int[i_lat]))

				varss_int /= weights
				# varss_int = f(latgrid)

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
						file.write(str(varss_int[i]*lat_facs[i]))
						file.write('\n')

					file.close()
			def create_manual_cloud_inputs():
				
				# Calculate clear sky fraction with random overlap assumption
				c1 = 0.
				c2 = 0.
				ctot = 0.
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


			interpdir = 'Input Data for RCM Interpolation/'
			outdir = project_dir+'Input Distributions/' #output file directory

			misr_cf_latalt_max = np.load(interpdir+'misr_cf_latalt.npy') # maximum-sized array of cloud data, to interpolate onto new smaller grid here
			misr_cf_latalt_max = np.nan_to_num(misr_cf_latalt_max)
			xbounds = np.linspace(-1,1,ncols+1)
			# latbounds = np.rad2deg(np.arcsin(xbounds))
			latbounds = np.linspace(-90,90,ncols+1)
			# latbounds = [-90,-66.5,-23.5,23.5,66.5,90]
			collats = np.zeros(ncols)
			colwidth = 180.0 / (ncols+1)

			for i in range(1,len(latbounds)):
				collats[i-1] = (latbounds[i-1] + latbounds[i])/2
				print(collats[i-1])

			latweights = np.cos(radians(collats))
			latgrid = collats
			latgridbounds = latbounds
			# pgrid = np.linspace(1000,1,nlays)
			# pgrid = np.linspace(1000,min_press,nlays)
			

			q_latp_max = np.load(interpdir+'q_latp.npy')
			o3_latp_max = np.load(interpdir+'o3_latp.npy')
			fal_lat_max = np.load(interpdir+'fal_lat.npy')
			t_latp_max = np.load(interpdir+'t_latp.npy')

			q_ps = np.load(interpdir+'q_ps.npy')
			o3_ps = np.load(interpdir+'o3_ps.npy')
			t_ps = np.load(interpdir+'t_ps.npy')

			q_lats = np.load(interpdir+'q_lats.npy')
			o3_lats = np.load(interpdir+'o3_lats.npy')
			fal_lats = np.load(interpdir+'fal_lats.npy')
			t_lats = np.load(interpdir+'t_lats.npy')

			pa = 0.3    
			sc = [1362.0]

			#pico2s = np.linspace(400e-6,3200e-6,num=5)

			#pico2s = np.logspace(-4,2,num=5,base=10.0)

			

			# pico2_facs = np.array([1e-2,0.0625,0.125,0.25,0.5,1,2,4,8])
			# if(forcing_expt==0):
			# 	gas_amt_fac_co2s = np.array([1,2])
			# else:
			# 	gas_amt_fac_co2s = np.array([2])
			# pico2_facs = np.array([1.0])
			pico2s = np.array([531e-6]) * gas_amt_fac_co2s
			# pico2 = 348e-6 #rcemip
			pin2s = [0.78084] #%
			pio2 = 0.20946 #%
			piar = 0.009340 * 0. #%
			pich4 = 1650e-9 * 0.  #rcemip

			# pico2s = np.array([420e-6])

			# pin2s = np.logspace(-1,2,num=10,base=10.0)
			# pin2s = [1.0]
			# pin2s = [2.0]

			# pico2s = [400e-6]


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
			sas = [0.2]
			#tboundms = np.linspace(250,335,num=10)
			tboundms = [288.4]

			#for pertcol in range(ncols):

			# global_lapses = np.linspace(-2,-8,5)
			global_lapses = [-5.7]

			cloud_loc_type = 0 # 0: pressure (hPa), 1: altitude (km), 2: temperature (K)    

			# cld_heights = np.linspace(0,12,5)
			cld_heights = [5.0]
			# cld_height = [5.0]
			# cld_taus = np.linspace(0.1,9.9,10)
			cld_taus = [9.9]

			# mixco2_prescribed_facs = np.array([0.03125,0.0625,0.125,0.25,0.5,1,2,4,8])
			mixco2_prescribed_fac = 1.0
			# mixco2_prescribed_facs = np.array([0.03125])
			psurf_override = 1000.
			# psurf_overrides = [1.]
			#fsws = np.linspace(200,500,num=8)

			# add_cld_alts = [0.0,6.1]
			add_cld_alts = [0.0]
			# lcs = np.linspace(10,2,10)
			# lcs = lcs * -1.

			# nperts = 5
			pert_thickness = 1000./nperts
			# pperts = np.linspace(1000,pert_thickness,nperts)
			# pperts = np.insert(pperts,0,np.array([2000.]),axis=0)
			# bottom up:
			pperts = np.linspace(1000-pert_thickness,0,nperts)
			# pperts = np.insert(pperts,0,np.array([2000.]),axis=0)
			# top down:
			# pperts = np.linspace(pert_thickness,1000,nperts)
			# pperts = np.insert(pperts,0,np.array([-1000.]),axis=0)
			# print pperts
			

			lf_as = [0.0] # 0.0 default
			# twarms = [288.,293.,298.,303.,308.]
			twarm = 288.
			# tcolds = [268.,263.,258.,253.,248.]
			tcold = 268.

			
			for mtransp_type in mtransp_types:
				for o3_source in o3_sources:
					i_lt = 0
					for lapse_type in lapse_types:
						i_h2osrc=0
						for h2o_source in h2o_sources:
							i_pco2=0
							for pico2 in pico2s:
								gas_amt_fac_co2 = gas_amt_fac_co2s[i_pco2]
								i_lfa = 0
								for lf_a in lf_as:
									lat_facs = 1.0 + abs(np.sin(np.deg2rad(collats))) * lf_a #multiply a variable by a latitude-dependent factor to change the meridional gradient
									lat_facs = lat_facs * ncols / sum(lat_facs)
									i_cf=0
									for gas_amt_fac_h2o in gas_amt_fac_h2os:
										i_lc = 0
										for lc in lcs:
											i_pp = 0
											for ppert in pperts:
												i_ch = 0
												for fsw in fsws:
													i_pso = 0
													for cld_tau in cld_taus:
														for sa in sas:
															for pin2 in pin2s:
																for gas_amt_fac_ch4 in gas_amt_fac_ch4s:
																	for gas_amt_fac_h2o in gas_amt_fac_h2os:

																		nloops = len(cld_heights)*len(fsws)\
																		*len(cld_taus)*len(tboundms)*len(sas)*len(pin2s)*len(pico2s)*len(lapse_types)*len(pperts)\
																		*ncols*nlays*ncloudcols*len(lf_as)*len(h2o_sources)*len(lcs)\
																		*len(gas_amt_fac_h2os)
																		print("Number of loops: ", nloops)
																		secsperloop = 0.5 #uni
																		# secsperloop = 1.5 #home
																		print("Estimated mins: ",nloops*secsperloop/60.)
																		print("pico2=",pico2)
																		print("h2o_source=",h2o_source)
																		print("lapse_type=",lapse_type)
																		print("o3_source=",o3_source)
																		print("mtransp_type",mtransp_type)


																		inv_sum = (pico2+pin2+pio2+piar+pich4) * 1e3

																		# nlays = nlayss[i_pso]
																		pgrid = np.linspace(psurf_override,min_press,nlays+1)
																		# pgrid = np.linspace(inv_sum,min_press,nlays+1)

																		add_cld_alt = add_cld_alts[i_pso]

																		#mnlcld
																		
																		
														
																		# sa = [sa] * ncols

																		if ( cloud_source == 0 ):
																			ncloudcols = shape(manual_clouds)[0]
																			create_manual_cloud_inputs()
																		elif ( cloud_source == 1 ):
																			create_misr_cloud_inputs(add_cld_alt)
																	
																		#pertvars = ['q','cc','ciwc','clwc','o3']
																		pertvars = ['q']
														
																		# shortnameslev = ['cc','clwc','o3','q','ciwc']
																		shortnameslev = ['q', 'o3','t']
																		longnameslev = {'cc':'Cloud fraction','clwc':'Cloud liquid water content (kg/kg)','o3':'Ozone mixing ratio','q':'Specific humidity (kg/kg)','ciwc':'Cloud ice water content (kg/kg)','t':'Temperature'}
																		#disttypelev = {'cc':'lat','clwc':'lat','o3':'lat','q':'lat','ciwc':'lat'}
																		#disttypelev = {'cc':'lat','clwc':'lat','o3':'lat','q':'lat','ciwc':'lat'}
																		disttypelev = {'cc':'lat','clwc':'lat','o3':'lat','q':'lat','ciwc':'lat','t':'lat'}
																
																	
																		shortnamessfc = ['fal']
																		longnamessfc = {'fal':'Surface albedo'}
																		disttypesfc = {'fal':'lat'}
																	
																		loop = 1
																	
																		print("Creating input files by interpolating ERA-Interim data")
																	
																		# interpolate_createprrtminput_lev('q',q_latp_max,q_ps,q_lats,lat_facs)
																		# interpolate_createprrtminput_lev('o3',o3_latp_max,o3_ps,o3_lats,[1.0]*ncols)
																		# interpolate_createprrtminput_sfc('fal',fal_lat_max,fal_lats,[1.0]*ncols)

																		interpolate_createprrtminput_lev('q',q_latp_max,q_ps,q_lats,[1.0]*ncols)
																		interpolate_createprrtminput_lev('o3',o3_latp_max,o3_ps,o3_lats,[1.0]*ncols)
																		# interpolate_createprrtminput_sfc('fal',fal_lat_max,fal_lats,[1.0]*ncols)
																		interpolate_createprrtminput_lev('t',t_latp_max,t_ps,t_lats,[1.0]*ncols)
																	
																		if (lapse_type == 2 or lapse_type==1):
																			lc = createlatdistbn('Doug Mason Lapse Rate vs Latitude')
																		
																		if(ncols>1 and lapse_type==0):
																			lc = [lc] * ncols
																		# print(lc, 'lc')
																		

																		lch = createlatdistbn('Cloud Top Height')
																		srh = createlatdistbn('Relative Humidity')
																		# srh = [0.8] * ncols
																		sa = createlatdistbn('Surface Reflectance')
																		# print(sa)
																		# sa = list(sa * lat_facs)
																		# sa = [0.2] * ncols
																		lcf = createlatdistbn('Cloud Fraction')
																		lcod = createlatdistbn('Cloud Optical Thickness')
																		tg = createlatdistbn('Surface Temperature')
																		# tg = createlatdistbn('Doug Mason Temperature vs Latitude')
																		# for i in range(len(tg)):
																		# 	tg[i] += 20.
																		# print(tg)
																		# tg = [290.] * ncols - abs(collats)
																		# tg = tg * lat_facs

																		# tg = [tboundm] * ncols
																		# tg = [270.594] * ncols

																	
																		#lc = [-5.88]*ncols
																		#lch = [4.46]*ncols
																		#srh = [0.787]*ncols
																		#lcf = [0.657]*ncols
																		#lcod = [3.978]*ncols
																	
																		mc = pico2  
																	
																		ur = ur_htr
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
																		# tp = 1.0
																		#fth = np.zeros(ncols)
																		#for i in range(ncols):
																		#    fth[i] = 15.0 - abs(collats[i])/18.0
																		fth = [500.] * ncols
																		ol = nlays
																		asp = 2.0   
																		cs = 0
																		pbo = 0 
																		fp = 0
																		ps1 = 0
																		af = 1.0																	
																		npb = 1
																		o3sw = 1
																		h2osw = 1
																		nl = nlays
																		asf = 4.0
																		tuf = 1.0   
																		htransp = 1.0 #reduce lapse rate to account for horizontal transport
																		ipe = 1
																		dp = 1
																		mtranspfac = 2.0
																		boxnetfluxfac = 0.2
																		# twarm = 288
																		# tcold = 268
																		phim = 45 * 3.14 / 180
																		ks = 0.25
																		kl = 0.25
																		eta = 0.75
																		planet_radius = 6.37e6
																		planet_rotation = 7.29e-5
																		t_min = 10.
																		sfc_heating = 0 #surface energy budget warms/cools surface? 1=yes, 0=no
																		playtype = 0 #pressure layer type. 0=equal p thickness, 1=sigma
																		
																		ur_toafnet = [4.0] * ncols
																		ur_seb = 1e10
																		couple_tgta = 1
																		mtranspon = 1
																		# gas_amt_fac_h2o = 1.0
																		# gas_amt_fac_co2 = 1.0
																		gas_amt_fac_o3 = 1.
																		# gas_amt_fac_ch4 = 1.0
																		# gas_amt_p_high_h2o = ppert
																		# gas_amt_p_low_h2o = ppert - pert_thickness
																		# gas_amt_p_high_h2o = 1e6
																		# gas_amt_p_low_h2o = 0.
																		gas_amt_p_high_h2o = 1000.
																		gas_amt_p_low_h2o = ppert #bottom up
																		# gas_amt_p_high_h2o = ppert
																		# gas_amt_p_low_h2o = 0.
																		gas_amt_p_high_co2 = 1e6
																		gas_amt_p_low_co2 = 0.
																		# gas_amt_p_high_co2 = ppert
																		# gas_amt_p_low_co2 = ppert - pert_thickness
																		# gas_amt_p_high_co2 = 1000.
																		# gas_amt_p_low_co2 = ppert
																		print(gas_amt_p_high_h2o,'to ', gas_amt_p_low_h2o)
																		gas_amt_p_high_o3 = 1e6
																		gas_amt_p_low_o3 = 0.
																		gas_amt_p_high_ch4 = 1e6
																		gas_amt_p_low_ch4 = 0
																		# gas_amt_p_high_ch4 = ppert
																		# gas_amt_p_low_ch4 = ppert-pert_thickness
																		gas_amt_pert_h2o = 1 #1 = on, 0=off
																		gas_amt_pert_co2 = 1 #1 = on, 0=off
																		gas_amt_pert_o3 = 1 #1 = on, 0=off
																		gas_amt_pert_ch4 = 1 #1 = on, 0=off
																		psurf_override = psurf_override #override the inventory base psurf calc and set explicit psurf; set < 0 to turn this option off.
																		mixco2_prescribed_on = 0
																		mixco2_prescribed = 400e-6 * mixco2_prescribed_fac
																		steps_before_toa_adj = 5
																		a_green = 0.4
																		b_green = 20.
																		c_green = 5.
																		H_green = 7.
																		cloudloctype = 2 #1 for altitude, 2 for pressure, 3 for temperature
																		surf_emiss_on = 1 #0 for no surface emission, 1 for normal surface emission
																		# lapse_type = 1
																		h2o_sb = 1 #h2o foreign broadening 0=off, 1=on
																		h2o_for = 1
																		# h2o_source = 2 # 0=ERA-I mixh2o, 1=MW67 RH, 2=Cess RH, 3=Kasting, 4=Ramirez
																		ur_mt = 1.0
																		# ur_mt = 4.
																		gas_addmolec_h2o = 0.0
																		# gas_addmolec_h2o = 1.0e19 * (ppert/1000.) * 3.83 / 3.3
																		# gas_addmolec_h2o = 3.82 / 2.31 * 1.0e19 * (ppert/1000.)**2
																		gas_addmolec_co2 = 0.0
																		gas_addmolec_o3 = 0.0
																		gas_addmolec_ch4 = 0.0
																		max_rh = 1.0 * 1e6
																	
																		ur1 = ur
																	
																		counter = 0
																	
																		ur = ur1
																	
																		params = [ncols,ncloudcols+1,pa,sc,list(tg),lc,days,mc,ur,icldm,rmin,hct,hcf,hcod,mct,mcf,mcod,lch,lcf,lcod,tp,sa,list(fth),ol,asp,cs,pbo,fswon,fsw,fp,srh,ps1,af,convecttype,
																		npb,o3sw,h2osw, nl, maxhtr, asf, tuf, pico2, pin2, pio2, htransp, ipe, dp, mtranspfac,boxnetfluxfac,pertlay,pertcol,list(collats),inversion_strength,inversion_col,
																		twarm,tcold,phim,ks,kl,eta,planet_radius,planet_rotation,list(latbounds),t_min,sebfac,sfc_heating,playtype,ur_htr,ur_toafnet,ur_seb,couple_tgta,mtranspon,min_press,
																		gas_amt_fac_h2o,gas_amt_fac_co2,gas_amt_fac_o3,gas_amt_p_high_h2o,gas_amt_p_low_h2o,gas_amt_p_high_co2,gas_amt_p_low_co2,gas_amt_p_high_o3,gas_amt_p_low_o3,
																		gas_amt_pert_h2o,gas_amt_pert_co2,gas_amt_pert_o3,psurf_override,mixco2_prescribed_on,mixco2_prescribed,steps_before_toa_adj,a_green,b_green,c_green,H_green,cloudloctype,
																		surf_emiss_on,lapse_type,h2o_for,h2o_sb,h2o_source,ur_mt,mtransp_type,steps_before_first_eqbcheck,gas_addmolec_h2o,gas_addmolec_co2,gas_addmolec_o3,max_rh,snapshot,
																		piar,pich4,gas_amt_p_high_ch4,gas_amt_p_low_ch4,gas_amt_pert_ch4,gas_amt_fac_ch4,gas_addmolec_ch4,forcing_expt,o3_source]
																		
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
																			#print('Loop # {}, {:4.1f}% complete'.format(counter,np.float(counter)/np.float(nloops)))
																	
																			if (p.returncode == 0):
																				break
																			elif (p.returncode == -11 or p.returncode == -8 or p.returncode == 11 or p.returncode == 12):
																				ur = ur*2
																			continue
																 
																		counter = counter + 1
												i_ch +=1
											i_pp+=1
										i_lc+=1
									i_cf+=1
								i_lfa+=1
							i_pco2+=1
						i_h2osrc+=1
					i_lt+=1

########################################################################################################################

os.system('say "Done"')
show()