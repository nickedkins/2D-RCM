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

plot_all_vert_profiles = 0
kluftfig=0
legends_on = 1
grids_on = 1

directories = [
'_Current Output/'
]

skip_ifn = []

directories = [
'/Users/nickedkins/Dropbox/GitHub Repositories/Uni/2D-RCM/_Useful Data/held82 replication/'
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
	plt.rcParams['font.size'] = 20
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

	plt.rcParams['lines.linewidth'] = 2.0 
	plt.rcParams['lines.markersize'] = 12

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
	# print "sinlats:", np.cos(np.deg2rad(lats))
	# print "sumlats:", sum(np.cos(np.deg2rad(lats)))
	x_avg = sum(x * np.cos(np.deg2rad(lats))) / sum(np.cos(np.deg2rad(lats)))
	print 'x_avg:', x_avg
	return x_avg

def latwghtavg_2d(x,lats):
	nlays = len(x[:,0])
	x_avg = np.zeros(nlays)
	for i in range(nlays):
		x_avg[i] = sum(x[i,:] * np.cos(np.deg2rad(lats))) / sum(np.cos(np.deg2rad(lats)))
	return x_avg

def latwghtavg_dim(x,lats,dim):
	x_avg = np.zeros(shape(x)[dim])
	print shape(x)
	for i in range(shape(x)[dim]):
		x_avg[i] = sum(x[i,:] * np.cos(np.deg2rad(lats))) / sum(np.cos(np.deg2rad(lats)))
	return x_avg

def make_csv(x,fname):
	f=open(str(fname+".csv"),"w+")
	for i in range(len(x)):
		mystr = str(str(i) + ',' + str(x[i]))
		f.write(mystr)
		f.write('\n')
	f.close
	


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


linestyles = ['-','--','-']

# colors = ['b','r','g','orange','purple','yellow','pink']



def readfile(fn,counter):

	output_file = directory + fn
	f = open(output_file,'r')
	params = f.readline().split()
	ncols=int(params[0])
	nlayersm=int(params[1])
	print nlayersm
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

tzm_store = np.zeros( ( 200, len(directories), 100 ) )
pico2_store = np.zeros( (len(directories), 100 ) )

cld_heights = np.linspace(1,10,10)


# pertfacs = [1.,-1.]

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
	ptrops_master = []
	ttrops_master = []
	ztrops_master = []
	tsurfs_master = []
	tavelm_master = []
	wklm2_master = []
	wklm3_master = []
	totdflumcols_master = []
	pavelm_master = []

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

	# nlayss = [25,50,100,200,400,600]

	# ncolss = [1,2,3,4,5,6,7,8,9,10]
	ncolss = [1,3,5,7,9]
	# ncolss = [1,2,4,6,8,10]

	fsws = [100,150,200,250,300,350]

	T1s = np.zeros(5)
	T2s = np.zeros(5)


	i_fn = 0
	for fn in a:
		if(i_fn in skip_ifn):
			i_fn +=1
			continue
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
		tavelm_master.append(tavelmcols)  
		wklm2_master.append(wklm2cols)
		wklm3_master.append(wklm3cols)
		totdflumcols_master.append(totdflumcols)
		pavelm_master.append(pavelmcols)  

		# T1s[int(i_fn)] = latwghtavg(tzmcols[0,:],boxlatcols[0,:])


		# if(mod(i_fn,2)==0):
		#     T1s[int(i_fn/2.)] = latwghtavg(tzmcols[0,:],boxlatcols[0,:])
		# if(mod(i_fn,2)==1):
		#     T2s[int(i_fn/2.)] = latwghtavg(tzmcols[0,:],boxlatcols[0,:])




		# plt.plot(ncolss[int(i_fn/2)],latwghtavg(tzmcols[0,:],boxlatcols[0,:]),'o',c='b')
		# plt.xlabel('Number of latitude columns')
		# plt.ylabel('Equilibrium $T_{surf}$ (K)')

		htrmlwcols = htrmcols[1:,:] - htro3cols - htrh2ocols

		conv_trop_ind_cols = np.zeros(ncols)

		for col in range(ncols):
			conv_trop_ind = int(convcols[i_fn,col])
			if conv_trop_ind > nlayersm:    
				conv_trop_ind = nlayersm
			conv_trop_ind_cols[col] = int(conv_trop_ind)

			p_trop = pzmcols[conv_trop_ind,col]
			t_trop = tzmcols[conv_trop_ind,col]
			z_trop = altzmcols[conv_trop_ind,col]
			tsurf = tzmcols[0,col]

			ptrops_master.append(p_trop)
			ttrops_master.append(t_trop)
			ztrops_master.append(z_trop)
			tsurfs_master.append(tsurf)


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

				conv_trop_ind = int(convcols[i_fn,col])
				if conv_trop_ind > nlayersm:    
					conv_trop_ind = nlayersm

				# print('DW LW at tropopause: ', totdflumcols[conv_trop_ind,col])
				# print('abs SW above tropopause: ', sum(A_oz_lcols[conv_trop_ind:nlayersm,col])*1362./4.)

				p_trop = pzmcols[conv_trop_ind,col]
				t_trop = tzmcols[conv_trop_ind,col]
				z_trop = altzmcols[conv_trop_ind,col]

				# print 'ptrop, ttrop, ztrop, psurf, tsurf'
				# print(p_trop,',', t_trop,',', z_trop/1000.,',', pzmcols[0,col],',', tzmcols[0,col])

				# plt.figure(2)
				# plt.plot(lapsecritcols[0,col],z_trop/1000.,'o')

				# for i in range(len(altzmcols[:,col])):
				#     print altzmcols[i,col]/1000., ',', pzmcols[i,col], ',', tzmcols[i,col]

				plt.figure(i1+1)
				# plt.figure(i_dir)
				plt.subplot(341)
				plt.title('tzm')
				# plt.plot(tzmcols[:,col],altzmcols[:,col]/1000.,ls=linestyles[i1],label='Spectral '+str(fn))
				# plt.semilogy(tzmcols[:,col],pzmcols[:,col],label=str(fn),ls=linestyles[i1])
				plt.semilogy(tavelmcols[:,col],pavelmcols[:,col],label=str(fn),ls=linestyles[i1])
				# plt.semilogy(tzmcols[:,col],pzmcols[:,col],'-o',label=str(fn))
				plt.xlabel('Temperature (K)')
				plt.ylabel('Pressure (hPa)')
				# plt.plot(tzmcols[:,col],altzmcols[:,col],'-o',label=str(fn))
				# plt.plot(tzmcols[conv_trop_ind,col],altzmcols[conv_trop_ind,col]/1000.,'*',markersize=20)
				plt.plot(tzmcols[conv_trop_ind,col],pzmcols[conv_trop_ind,col],'o',alpha=0.5)
				plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))
				if(grids_on==2):
					plt.gca().minorticks_on()
				if(grids_on==1):
					plt.grid(which='both',axis='both')
				if(legends_on==1):
					plt.legend()    
				
				plt.figure(i1+1)
				plt.subplot(342)
				plt.title('htrm')
				# plt.semilogy(htrmcols[:,col],pzmcols[:,col],ls=linestyles[i1],label=str(fn))
				plt.semilogy(htrmcols[:nlayersm-1,col],pzmcols[:nlayersm-1,col],'-o',label=str(fn))
				# plt.xlim(-5,5)
				plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))
				plt.axvline(-0.03,ls='--')
				plt.axvline(0.03,ls='--')
				if(legends_on==1):
					plt.legend()
				
				plt.figure(i1+1)
				plt.subplot(343)
				plt.title('totuflum')
				plt.semilogy(totuflumcols[:,col],pzmcols[:,col],ls=linestyles[i1],label='up '+dir_label)
				plt.semilogy(totdflumcols[:,col],pzmcols[:,col],ls=linestyles[i1],label='down '+dir_label)
				plt.semilogy(totdflumcols[:,col]-totuflumcols[:,col],pzmcols[:,col],ls=linestyles[i1],label='net '+dir_label)
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
				# plt.xlim(-5,5)
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
				# plt.xlim(-5,5)
				plt.ylim(max(pzmcols[:,col]),min(pzmcols[:,col]))
				plt.axvline(-0.03,ls='--')
				plt.axvline(0.03,ls='--')
				if(legends_on==1):
					plt.legend()

				i3+=1

		i_fn += 1


	# plt.plot(ncolss,T2s-T1s,'-o',mec='none')
	# # plt.axhline(1.2,linestyle='--')
	# plt.xlabel('Number of latitude columns')
	# plt.ylabel('Change in temperature from double CO$_2$')

	# plt.plot(ncolss,T1s+14.,'-o')
	# # plt.plot(ncolss,T2s,'-o')
	# plt.xlabel('Number of latitude columns')
	# plt.ylabel('Temperature (K)')


	# i1 += 1
	i_dir+=1

	# plt.legend()

	# plt.plot(ncolss,T1s)
	


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
	ptrops_master=np.array(ptrops_master)
	ttrops_master=np.array(ttrops_master)
	ztrops_master=np.array(ztrops_master)
	tsurfs_master=np.array(tsurfs_master)
	tavelm_master = np.array(tavelm_master)
	wklm2_master=np.array(wklm2_master)
	wklm3_master=np.array(wklm3_master)
	totdflumcols_master=np.array(totdflumcols_master)
	pavelm_master = np.array(pavelm_master)
	
	# master indices for conv_trop_ind: master[file][column]
	# to get trop values for tzm, use tzm_master[file,conv_trop_ind_master[file,column],column][0][file] (slightly ugly but functional)
	conv_trop_ind_master = np.array(conv_trop_ind_master)

	lats = boxlatcols_master[0,0,:]

	ttrops = np.zeros((6,ncols))

	# make_csv(tzm_master[0,:,0],'tzm')
	# make_csv(tavelm_master[0,:,0],'tavelm')
	# make_csv(wklm1_master[0,:,0],'wklm1 (h2o)')
	# make_csv(wklm2_master[0,:,0],'wklm2 (co2)')
	# make_csv(wklm3_master[0,:,0],'wklm3 (o3)')
	# make_csv(wbrodlm_master[0,:,0],'wbrodlm')
	# make_csv(totuflumcols_master[0,:,0],'totuflum')
	# make_csv(totdflumcols_master[0,:,0],'totdflum')
	# make_csv(pavelm_master[0,:,0],'pavelm')
	# make_csv(pzm_master[0,:,0],'pzm')
	# make_csv(altzm_master[0,:,0],'altzm')
	# make_csv(abspncols_master[0,:,0]*1362./4.,'abs sw h2o')
	# make_csv(A_oz_lcols_master[0,:,0]*1362./4.,'abs sw o3')

	# print sum(abspncols_master[0,:,0]) * 413.177 + sum(A_oz_lcols_master[0,:,0]) * 413.177 + abs_surf_cols_master[0,0,0]
	# print abs_surf_cols_master[0,0,0]
	# print sum(A_oz_lcols_master[0,:,0]) * 413.177
	# print sum(abspncols_master[0,:,0]) * 413.177


	lcs = np.arange(-10,-3)

	data=np.genfromtxt('/Users/nickedkins/Dropbox/Figure Data/held82_x[lapse]_y[ztrop].txt',delimiter=',')
	lapse_h82=data[:,0]*-1.
	ztrop_h82=data[:,1]

	plt.figure(1)
	plt.subplot(221)
	plt.plot(lcs,ptrops_master,'-o')
	plt.xlabel('Lapse Rate (K/km)')
	plt.ylabel('Tropopause Pressure (hPa)')
	plt.subplot(222)
	plt.plot(lcs,ttrops_master,'-o')
	plt.xlabel('Lapse Rate (K/km)')
	plt.ylabel('Tropopause Temperature (K)')
	plt.subplot(223)
	plt.plot(lcs,tzm_master[:,0,:],'-o')
	plt.xlabel('Lapse Rate (K/km)')
	plt.ylabel('Surface Temperature (K)')
	plt.subplot(224)
	plt.plot(lcs,ztrops_master*1e-3,'-o')
	plt.plot(lapse_h82,ztrop_h82,'-',label='Held 1982')
	plt.legend()
	plt.xlabel('Lapse Rate (K/km)')
	plt.ylabel('Tropopause Height (km)')


	# i_co2 = 0
	# for i_co2 in range(6):
	# 	i_lat = 0
	# 	for i_lat in range(ncols):
	# 		cti = conv_trop_ind_master[i_co2,i_lat]
	# 		ttrops[i_co2,i_lat] = tzm_master[i_co2,cti,i_lat]




	# print np.diagonal(tzm_master[:,conv_trop_ind_master[:,0],0])
	# print(pzm_master[:,conv_trop_ind_master[0,:],0])
	
	# _kluftfig
	if(kluftfig==1):
		plt.figure(i_dir)
		plt.gcf().suptitle(dir_label)

		# if(i_fn==1):

		lats = boxlatcols_master[0,0,:]

		for col in range(ncols):

			plt.subplot(121)
			plt.plot(tzm_master[:,0,col]-tzm_master[1,0,col],np.diagonal(pzm_master[:,conv_trop_ind_master[:,0],col])-np.diagonal(pzm_master[:,conv_trop_ind_master[:,0],col])[1],'-',alpha=0.5,color='grey')
			plt.xlabel("$\Delta T_{surf}$",labelpad=20)
			plt.ylabel("$\Delta p_{trop}$")
			plt.legend()
			plt.grid(True)

			plt.subplot(122)
			# plt.plot(tzm_master[:,0,col]-tzm_master[1,0,col],np.diagonal(tzm_master[:,conv_trop_ind_master[:,0],col])-np.diagonal(tzm_master[:,conv_trop_ind_master[:,0],col])[1],'-',alpha=0.5,color='grey')
			# plt.plot(tzm_master[:,0,col]-tzm_master[1,0,col],ttrops[:,col]-ttrops[1,col],'-',alpha=0.5,color='grey')
			plt.plot(tzm_master[:,0,col]-tzm_master[1,0,col],ttrops[:,col]-ttrops[1,col],'-',alpha=0.5,label=col)
			print tzm_master[:,0,col]-tzm_master[1,0,col]

		# for i_fn in range(len(a)):
		# 	plt.plot(latwghtavg(tzm_master[i_fn,0,:]-tzm_master[1,0,:],lats ),latwghtavg( ttrops[i_fn,:]-ttrops[1,:],lats), 'o',c='b' )



		plt.subplot(121)
		
		# plt.plot(latwghtavg_2d(tzm_master[:,0,:],lats)-latwghtavg_2d(tzm_master[1,0,:],lats),latwghtavg_2d(np.diagonal(pzm_master[:,conv_trop_ind_master[:,0],:]), lats ) )

		plt.subplot(122)
		data = np.genfromtxt('/Users/nickedkins/Dropbox/Figure Data/2020-01 (January)/x[dTsurf]_y[dTtrop]_kluft2019fig8_mw67.txt',delimiter=',')
		x_k19 = data[:,0]
		y_k19 = data[:,1]
		plt.plot(x_k19,y_k19,'o',label='Kluft using Manabe 67')
		data = np.genfromtxt('/Users/nickedkins/Dropbox/Figure Data/2020-01 (January)/x[dTsurf]_y[dTtrop]_kluft2019fig8_phat.txt',delimiter=',')
		x_k19 = data[:,0]
		y_k19 = data[:,1]
		plt.plot(x_k19,y_k19,'--',c='b',alpha=0.5,label='PHAT')
		plt.plot(tzm_master[:,0,0]-tzm_master[1,0,0],tzm_master[:,0,0]-tzm_master[1,0,0],ls='--',c='r',alpha=0.5,label='FAP')       
		plt.plot(latwghtavg_dim(tzm_master[:,0,:]-tzm_master[1,0,:],lats,1), latwghtavg_dim(ttrops[:,:]-ttrops[1,:],lats,1)  ,'-o',label='2D Average')
		plt.xlabel("$\Delta T_{surf}$",labelpad=20)
		plt.ylabel("$\Delta T_{trop}$")
		plt.grid(True)
		plt.legend()

	# else:

	# 	lats = boxlatcols_master[0,0,:]

	# 	for i_co2 in range(6):

	# 		plt.subplot(122)
	# 		# plt.plot(latwghtavg_dim(tzm_master[:,0,:]-tzm_master[1,0,:],lats,1), latwghtavg_dim(ttrops[:,:]-ttrops[1,:],lats,1)  ,'-o',label='1D Average')
	# 		plt.plot(latwghtavg(tzm_master[i_co2,0,:]-tzm_master[1,0,:],lats), latwghtavg( ttrops[i_co2,:] - ttrops[1,:],lats ), 'o',c='red')
	# 		plt.xlabel("$\Delta T_{surf}$",labelpad=20)
	# 		plt.ylabel("$\Delta T_{trop}$")
	# 		plt.grid(True)
	# 		plt.legend()



	box_abssw_tot_master = abs_surf_cols_master + abs_h2o_cols_master + abs_o3_cols_master
	# boxtotnetflux_master = meridtransp_master[:,0,:] + box_abssw_tot_master[:,0,:] - totuflumcols_master[:,-1,:]

	

	# Bs = pzm_master[:,conv_trop_ind_master[:,range(ncols)],range(ncols)]
	As = np.zeros(len(a))
	Bs = np.zeros(len(a))

	# for i_fn in range(len(a)):
	#     As[i_fn] = tzm_master[i_fn,conv_trop_ind_master[i_fn],:] / tzm_master[0,conv_trop_ind_master[0],:]
	#     Bs[i_fn] = ( pzm_master[i_fn,conv_trop_ind_master[i_fn],:] / pzm_master[0,conv_trop_ind_master[0],:] ) ** (lapsecritcols_master[i_fn,0,0]*0.287/9.8)


	pico2_facs = np.array([1./32.,1./16.,1./4.,1,4,16,32])
	pico2s = np.array([420e-6]) * pico2_facs

	# plt.title('Absolute value of influence of $p_{trop}$ and $T_{trop}$ on $T_{surf}$')
	# plt.title('Influence of $p_{trop}$ and $T_{trop}$ on $T_{surf}$')
	# # plt.plot(pico2s,abs((As-1.0)*tzm_master[0,0,:]),'-o',label='$T_{trop}$')
	# # plt.plot(pico2s,abs((Bs-1.0)*tzm_master[0,0,:]),'-o',label='$p_{trop}$')
	# plt.semilogx(pico2s,(As-1.0)*tzm_master[0,0,:],'-o',label='$T_{trop}$')
	# plt.plot(pico2s,(Bs-1.0)*tzm_master[0,0,:],'-o',label='$p_{trop}$')
	# # plt.plot(pico2s,(As*Bs-1.0)*tzm_master[0,0,:],'-o',label='AB (trop)')
	# plt.axhline(0,linestyle='--')
	# plt.xlabel('$p_{I,CO_2}$ (bar)')
	# plt.ylabel('Change in $T_{surf}$ caused by change in variable (K)')
	# plt.legend()

	# pperts = np.linspace(1000,50,20)
	# pperts = np.insert(pperts,0,np.array([2000.]),axis=0)

	# qfacs = [0.33,0.5,1.0,2.0,3.0]
	# qfacs = np.linspace(0.25,4.0,10)
	# lapses = np.linspace(10,2,10)
	# xlims = [lapses[0],lapses[-1]]
	# # ylims = [tzm_master[0,0,0]-tzm_master[3,0,0],tzm_master[-1,0,0]-tzm_master[3,0,0]]
	# ylims = np.array([tzm_master[0,0,0],tzm_master[-1,0,0]])
	# # plt.plot(qfacs,tzm_master[:,0,0]-tzm_master[3,0,0],'-o')
	# plt.plot(lapses,tzm_master[:,0,0]+6.,'-o')
	# plt.plot(xlims,ylims+6.,'--')
	# plt.xlabel('Lapse rate (K/km)',labelpad=20)
	# plt.ylabel('Temperature (K)',labelpad=20)

	nperts = 10
	pert_thickness = 1000./nperts
	pperts = np.linspace(1000,pert_thickness,nperts)

	# pperts = np.linspace(1000-pert_thickness,0,nperts)

	# pperts = np.insert(pperts,0,np.array([2000.]),axis=0)
	
	# bottom up:
	# if(directions[i_dir-1]==1):
	#     pperts = np.linspace(1000,pert_thickness,nperts)
	#     # pperts = np.insert(pperts,0,np.array([2000.]),axis=0)
	# # top down:
	# if(directions[i_dir-1]==2):
	#     pperts = np.linspace(pert_thickness,1000,nperts)
	#     # pperts = np.insert(pperts,0,np.array([0.]),axis=0)


	# print 'Total change: ', sum(tzm_master[:,0,0]-tzm_master[0,0,0]), 'K'

	ncloudcolss = [1,2,3,4,5,6,7]
	lcs = [-10,-8,-6,-4,-2,0]
	qs = [0.25,0.5,1.,2.,4.,8.]
	sas = [0.,0.1,0.2,0.3,0.4]

	# plt.plot(sas,tzm_master[:,0,:],'-o')
	# plt.xlabel('sa')
	# plt.ylabel('Equilibrium $T_{surf}$ (K)')
	

	# plt.plot(ncloudcolss,tzm_master[:,0,:],'-o')
	# plt.xlabel('Number of independent cloud columns')
	# plt.ylabel('Equilibrium $T_{surf}$ (K)')

	# print len(tzm_master[1:,0,0]-tzm_master[0,0,0]),len(pperts)

	# print tzm_master[:,0,0]

	# plt.figure(1)
	# plt.title('Change in equilibrium surface temperature with an \n increase in number of H$_2$O molecules in the 100 hPa above a given pressure')
	# plt.plot((tzm_master[1:,0,0]-tzm_master[0,0,0]),pperts,'-o',label=dir_label)
	# plt.xlabel('$\Delta T$ (K)')
	# plt.ylabel('Pressure at bottom of perturbation (hPa)')
	# plt.ylim(1000,0)
	# plt.legend()
	# # plt.xlim(0,0.2)

	# for i in range(1,len(a)):
	#     print pperts[i-1],',', tzm_master[i,0,0]-tzm_master[0,0,0]

	# pperts = [2000.,1000.]

	# pperts = [2000.,1000.,800.,600.,400.,350.,300.,250.,200.,150.,100.]

	# print 'dT = ', latwghtavg(tzm_master[1,0,:] - tzm_master[0,0,:],boxlatcols_master[0,0,:])
	# print 'T1 = ',latwghtavg(tzm_master[0,0,:],boxlatcols_master[0,0,:])
	# print 'T2 = ',latwghtavg(tzm_master[1,0,:],boxlatcols_master[0,0,:])

	# plt.figure(1)
	# # dT = latwghtavg(tzm_master[1,0,:] - tzm_master[0,0,:],boxlatcols_master[0,0,:])
	# # plt.plot(ncolss[i_dir-1],dT,'o')
	# plt.plot(boxlatcols_master[0,0,:], tzm_master[1,0,:]-tzm_master[0,0,:],'-o',label='Columns: '+str(ncolss[i_dir-1]))
	# plt.xlabel('Latitude')
	# plt.ylabel('$\Delta T_g$')
	# plt.legend()
	# plt.ylim(0)


	# tzm_master = tzm_master - 20.
	# y_endpoints=[tzm_master[0,0,0],tzm_master[-1,0,0]]
	# x_endpoints = [lapsecritcols_master[0,0,0],lapsecritcols_master[-1,0,0]]
	# plt.figure(1)
	# plt.plot(lapsecritcols_master[:,0,0],tzm_master[:,0,0],'-o')
	# plt.plot(x_endpoints,y_endpoints,'--')
	# plt.xlabel('Critical lapse rate (K/km)')
	# plt.ylabel('Surface temperature (K)')


	# tzm_master = tzm_master - 20.
	# h2o_facs = np.linspace(1.,10.,10)
	# y_endpoints=[tzm_master[0,0,0],tzm_master[-1,0,0]]
	# x_endpoints = [h2o_facs[0],h2o_facs[-1]]
	# plt.figure(1)
	# plt.plot(h2o_facs,tzm_master[:,0,0],'-o')
	# plt.plot(x_endpoints,y_endpoints,'--')
	# plt.xlabel('H$_2$O factor')
	# plt.ylabel('Surface temperature (K)')

	# tzm_master = tzm_master - 20.
	# cld_taus = np.linspace(0,9.9,10)
	# y_endpoints=[tzm_master[1,0,0],tzm_master[-1,0,0]]
	# x_endpoints = [cld_taus[1],cld_taus[-1]]
	# plt.figure(1)
	# plt.plot(cld_taus[1:],tzm_master[1:,0,0],'-o')
	# plt.plot(x_endpoints,y_endpoints,'--')
	# plt.xlabel('Cloud $\tau$')
	# plt.ylabel('Surface temperature (K)')
	# plt.figure(2)
	# plt.plot(cld_taus[1:],totuflumcols_master[1:,-1,0])
	# plt.plot(cld_taus[1:],box_abssw_tot_master[1:,-1,0])

	# init_h2o_molec = sum(wklm1_master[0,:,0])

	# print('Integrated temperature change: {:4.2f} K'.format(sum(tzm_master[:,0,0]-tzm_master[0,0,0])))
	# print('Initial total H2O molecules: {:4.2e}'.format(init_h2o_molec))
	# print('Integrated water vapor change (total molecules): {:4.2e}'.format(sum(wklm1_master[1:,:,0]-wklm1_master[0,:,0])))
	# print('Integrated water vapor change (additional molecs as % of original): {:0.0f}%'.format( sum( wklm1_master[:,:,0]-wklm1_master[0,:,0] )/init_h2o_molec*100. ))
	# for i in range(len(a)):
	#     print i, sum(wklm1_master[i,:,0])/init_h2o_molec
	# print

	


	# for i_fn in range(len(a)):

	#     for col in range(ncols):
	#         plt.figure(1)
	#         plt.subplot(121)
	#         plt.plot(rel_hum_master[i_fn,:,col]*100.,pzm_master[i_fn,1:,0],'-o',label='lat: '+str(int(boxlatcols_master[i_fn,0,col])))
	#         plt.xlabel('Relative Humidity (%)')
	#         plt.ylabel('Pressure (hPa)')
	#         plt.ylim(1000,0)
	#         plt.axvline(100,linestyle='--')
	#         plt.legend()

	#         plt.subplot(122)
	#         plt.plot(tzm_master[i_fn,:,col],pzm_master[i_fn,:,0],'-o',label='lat: '+str(int(boxlatcols_master[i_fn,0,col])))
	#         plt.ylim(1000,0)
	#         plt.xlabel('Temperature (K)')
	#         plt.ylabel('Pressure (hPa)')
	#         plt.legend()

	# for i_fn in range(len(a)):

	#     # vabsmax = np.max(abs(tzm_master[0,:,:]-tzm_master[1,:,:]))

	#     # plt.figure(1)

	#     # plt.subplot(211)
	#     # plt.contourf(wklm1_master[i_fn,:,:],20,vmax=6e21)
	#     # plt.colorbar()

	#     plt.figure(i_fn)
	#     plt.subplot(121)
	#     plt.contourf(wklm1_master[i_fn,:,:],10)
	#     plt.colorbar()
	#     plt.subplot(122)
	#     plt.contourf(rel_hum_master[i_fn,:,:],10,vmax=1)
	#     plt.colorbar()
		
		# plt.contourf(tzm_master[0,:,:]-tzm_master[1,:,:],cmap='bwr',vmax=vabsmax,vmin=-vabsmax)
		

	myvmin=0
	myvmax=10
	mymp=1


	# titles = ['ERA-Interim','Manabe-Wetherald','Cess','Kasting-Ackerman','Ramirez']

	# for i_fn in [0,2,4,6,8]:

	#     print i_fn

		# XX,YY = np.meshgrid(boxlatcols_master[i_fn,0,:],pzm_master[i_fn,1:,0])

		# plt.figure(i_fn)
		# plt.title(titles[i_fn])
		# plt.contourf(XX,YY,wklm1_master[i_fn,:,:]/wklm1_master[0,:,:],100,cmap='bwr',norm=MidpointNormalize(midpoint=mymp,vmin=myvmin, vmax=myvmax))
		# plt.gca().set_yscale('log')
		# plt.ylim(1000,1000/60.)
		# plt.xlabel('Latitude')
		# plt.ylabel('Pressure')
		# cb = plt.colorbar()
		# cb.set_label(r'Number of H$_2$O molecules',rotation=270,labelpad=20)

		# plt.figure(1)
		# # print latwghtavg_2d(wklm1_master[i_fn,:,:],boxlatcols_master[i_fn,0,:])
		# plt.loglog(latwghtavg_2d(wklm1_master[i_fn,:,:]/wbrodlm_master[i_fn,:,:],boxlatcols_master[i_fn,0,:]),pzm_master[i_fn,1:,0],'-o',label=titles[i_fn])
		# plt.ylim(1000,50)
		# plt.xlabel('H$_2$O mixing ratio')
		# plt.ylabel('Pressure (hPa)')
		# plt.xlim(3e-6)
		# # plt.gca().minorticks_on()
		# plt.grid(which='both',axis='both')


		# plt.figure(1)
		# plt.semilogy(latwghtavg_2d(wklm1_master[i_fn+1,:,:]/wbrodlm_master[i_fn+1,:,:]-wklm1_master[i_fn,:,:]/wbrodlm_master[i_fn,:,:],boxlatcols_master[i_fn,0,:]),pzm_master[i_fn,1:,0],'-o',label=titles[i_fn/2])
		# plt.ylim(1000,50)
		# plt.xlabel('Change in H$_2$O mixing ratio for double CO$_2$')
		# plt.ylabel('Pressure (hPa)')
		# # plt.xlim(3e-6)
		# plt.grid(which='both',axis='both')

		# plt.figure(1)
		# plt.minorticks_on()
		# plt.grid(which='both',axis='both')
		# plt.semilogy(latwghtavg_2d(rel_hum_master[i_fn+1,:,:]-rel_hum_master[i_fn,:,:],boxlatcols_master[i_fn,0,:]),pzm_master[i_fn,1:,0],'-o',label=titles[i_fn/2])
		# plt.ylim(1000,50)
		# plt.xlabel('Change in relative humidity for double CO$_2$')
		# plt.ylabel('Pressure (hPa)')

		# plt.figure(1)
		# # plt.minorticks_on()
		# plt.grid(which='both',axis='both')
		# plt.semilogy(latwghtavg_2d(tzm_master[i_fn+1,:,:]-tzm_master[i_fn,:,:],boxlatcols_master[i_fn,0,:]),pzm_master[i_fn,:,0],'-o',label=titles[i_fn/2])
		# plt.ylim(1000,50)
		# plt.xlim(-3,3)
		# plt.axvline(0,linestyle='--')
		# plt.xlabel('Change in temperature for double CO$_2$ (K)')
		# plt.ylabel('Pressure (hPa)')

		# plt.figure(1)
		# # plt.minorticks_on()
		# plt.grid(which='both',axis='both')
		# plt.semilogy(latwghtavg_2d(tzm_master[i_fn,:,:],boxlatcols_master[i_fn,0,:]),pzm_master[i_fn,:,0],'-o',label=titles[i_fn/2])
		# plt.ylim(1000,50)
		# # plt.xlim(-3,3)
		# plt.axvline(0,linestyle='--')
		# plt.xlabel('Temperature (K)')
		# plt.ylabel('Pressure (hPa)')

		# plt.figure(1)
		# # print latwghtavg_2d(wklm1_master[i_fn,:,:],boxlatcols_master[i_fn,0,:])
		# plt.semilogy(latwghtavg_2d(tzm_master[i_fn,:,:],boxlatcols_master[i_fn,0,:]),pzm_master[i_fn,:,0],'-o',label=titles[i_fn])
		# plt.ylim(1000,50)
		# plt.xlabel('Temperature (K)')
		# plt.ylabel('Pressure (hPa)')
		# # plt.xlim(3e-6)
		# # plt.gca().minorticks_on()
		# plt.grid(which='both',axis='both')




	# plt.figure(1)
	# plt.title('ERA-Interim')
	# plt.contourf(tzm_master[1,:,:]-tzm_master[0,:,:],20,cmap='bwr',norm=MidpointNormalize(midpoint=mymp,vmin=myvmin, vmax=myvmax))
	# plt.colorbar()

	# # divnorm = colors.DivergingNorm(vmin=-500, vcenter=0, vmax=4000)

	# plt.figure(2)
	# plt.title('Manabe-Wetherald')
	# plt.contourf(tzm_master[3,:,:]-tzm_master[2,:,:] - (tzm_master[1,:,:]-tzm_master[0,:,:]),20,cmap='bwr',norm=MidpointNormalize(midpoint=mymp,vmin=myvmin, vmax=myvmax))
	# plt.colorbar()

	# plt.figure(3)
	# plt.title('Cess')
	# plt.contourf(tzm_master[5,:,:]-tzm_master[4,:,:] - (tzm_master[1,:,:]-tzm_master[0,:,:]),20,cmap='bwr',norm=MidpointNormalize(midpoint=mymp,vmin=myvmin, vmax=myvmax))
	# plt.colorbar()

	# plt.figure(4)
	# plt.title('Kasting-Ackerman')
	# plt.contourf(tzm_master[7,:,:]-tzm_master[6,:,:]- (tzm_master[1,:,:]-tzm_master[0,:,:]),20,cmap='bwr',norm=MidpointNormalize(midpoint=mymp,vmin=myvmin, vmax=myvmax))
	# plt.colorbar()

	# plt.figure(5)
	# plt.title('Ramirez')
	# plt.contourf(tzm_master[9,:,:]-tzm_master[8,:,:]- (tzm_master[1,:,:]-tzm_master[0,:,:]),20,cmap='bwr',norm=MidpointNormalize(midpoint=mymp,vmin=myvmin, vmax=myvmax))
	# plt.colorbar()

	# i_fn = 0
	# for fn in a:
	#     # if (fn=='.DS_Store'):
	#     #     continue
	#     plt.figure(1)
	#     plt.subplot(221)
	#     plt.plot(boxlatcols_master[i_fn,0,:],meridtransp_master[i_fn,0,:],'-o',label=str(fn))
	#     plt.axhline(0)
	#     plt.xlabel('Latitude')
	#     plt.ylabel('mtransp')
	#     plt.legend()

	#     plt.subplot(222)
	#     plt.plot(boxlatcols_master[i_fn,0,:],tzm_master[i_fn,0,:],'-o',label=str(fn))
	#     # plt.axhline(0)
	#     plt.xlabel('Latitude')
	#     plt.ylabel('temp')
	#     plt.legend()

	#     plt.figure(1)
	#     plt.subplot(223)
	#     plt.plot(boxlatcols_master[i_fn,0,:],boxtotnetflux_master[i_fn,:],'-o',label=str(fn))
	#     plt.axhline(0)
	#     plt.xlabel('Latitude')
	#     plt.ylabel('totflux')
	#     plt.legend()

	#     i_fn+=1




	# filenames = np.array(filenames[0])

	# x = np.linspace(-90,90,ncols)

	# lats = boxlatcols_master[:,0,:]

	abs_sw_tot = abs_h2o_cols_master + abs_o3_cols_master + abs_surf_cols_master

	# for i_fn in range(0,len(boxlatcols_master[:,0,0])):

	#     plt.plot(boxlatcols_master[i_fn,0,:],meridtransp_master[i_fn,0,:],'-o')
	#     # plt.plot(boxlatcols_master[i_fn,0,:],abs_sw_tot[i_fn,0,:] - totuflumcols_master[i_fn,-1,:],'-o')
	#     plt.axhline(0)

	# twarms = [288.,293.,298.,303.,308.]
	tcolds = [268.,263.,258.,253.,248.]

	# for i_fn in range(len(a)):

	#     # plt.subplot(222)
	#     plt.plot(boxlatcols_master[i_fn,0,:],meridtransp_master[i_fn,0,:],'-o')
	#     plt.axhline(0)
	#     plt.xlabel('Latitude')
	#     plt.ylabel('mtransp')



	# print tzm_master[:,0,:], boxlatcols_master[0,0,:], latwghtavg(tzm_master[:,0,:],boxlatcols_master[0,0,:])
	
	# for i_fn in range(len(a)):
	#     fig=plt.figure(1)
		
	#     plt.subplot(231)
	#     plt.plot(boxlatcols_master[i_fn,0,:],tzm_master[i_fn,0,:]-tzm_master[0,0,:],'-o',label=a[i_fn])
	#     plt.xlabel('Latitude')
	#     plt.ylabel('$\Delta T$')
	#     plt.legend()

	#     plt.subplot(232)
	#     plt.plot(boxlatcols_master[i_fn,0,:],meridtransp_master[i_fn,0,:]-meridtransp_master[0,0,:],'-o',label=a[i_fn])
	#     plt.xlabel('Latitude')
	#     plt.ylabel('$\Delta$ meridional transport')
	#     plt.legend()

	#     plt.subplot(233)
	#     plt.plot(boxlatcols_master[i_fn,0,:],pzm_master[i_fn,conv_trop_ind_master[i_fn,range(ncols)],range(ncols)]-pzm_master[0,conv_trop_ind_master[0,range(ncols)],range(ncols)],'-o',label=a[i_fn])
	#     plt.xlabel('Latitude')
	#     plt.ylabel('$\Delta $Tropopause pressure (hPa)')
	#     plt.legend()

	#     plt.subplot(234)
	#     plt.plot(boxlatcols_master[i_fn,0,:],tzm_master[i_fn,conv_trop_ind_master[i_fn,range(ncols)],range(ncols)]-tzm_master[0,conv_trop_ind_master[0,range(ncols)],range(ncols)],'-o',label=a[i_fn])
	#     plt.xlabel('Latitude')
	#     plt.ylabel('$\Delta $Tropopause temperature (K)')
	#     plt.legend()

	#     plt.subplot(235)
	#     plt.plot(boxlatcols_master[i_fn,0,:],lapsecritcols_master[i_fn,0,:]-lapsecritcols_master[0,0,:],'-o',label=a[i_fn])
	#     plt.xlabel('Latitude')
	#     plt.ylabel('$\Delta$ Lapse rate (K/km)')
	#     plt.legend()

		# fig.suptitle('Increased surface albedo in tropics, decreased at poles',y=1.0)


	# plt.figure(1)
	# plt.subplot(221)
	# plt.plot(lats,tzm_master[0,0,range(ncols)])
	# plt.plot(lats,tzm_master[1,0,range(ncols)])
	# plt.subplot(222)  
	# plt.plot(lats,tzm_master[0,conv_trop_ind_master[0,:],range(ncols)])
	# plt.plot(lats,tzm_master[1,conv_trop_ind_master[1,:],range(ncols)])
	# plt.legend()


	# if (dir_label == 'Manabe-Wetherald 3 Clouds' or dir_label == 'Manabe-Wetherald Warmer'):
	#     pperts = np.linspace(1000,0,10)

	# plt.figure(1)
	# for i_fn in range(len(a)-1):
	#     if(a[i_fn]=='.DS_Store'):
	#         continue

	# plt.figure(1)
	# plt.title('Change in surface temperature with a perturbation \n in H$_2$O mixing ratio in each 50 hPa range')
	# plt.plot((tzm_master[1:,0,:]-tzm_master[0,0,:]),pperts,'-o',label=dir_label)
	# plt.xlabel('$\Delta T_{surf}$ (K)')
	# plt.ylabel('Pressure at bottom of 50 hPa H$_2$O perturbation region (hPa)')
	# plt.ylim(1000,0)

	# print 'Total delta T: {:4.2f}'.format(sum(tzm_master[1:,0,:]-tzm_master[0,0,:]))


# for i in range( shape(tzm_master)[0] ):

#     lats = boxlatcols_master[i,0,:]
#     pzms = pzm_master[i,:,0]
#     tzms = tzm_master[i,:,:]
#     altzms = altzm_master[i,:,0]
	# tsgm_weighted = sum(tzms[0,:] * np.cos(np.deg2rad(lats)) / sum(np.cos(np.deg2rad(lats)) ))



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
print('Done')
plt.tight_layout()
show()