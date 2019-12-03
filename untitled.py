import numpy as np
# import matplotlib as matplotlib
import matplotlib.pyplot as plt
from pylab import *
from scipy import interpolate
from os import listdir
import matplotlib.colors as colors
# import pandas as pd

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

data1 = np.genfromtxt('/Users/nickedkins/Dropbox/Figure Data/x[dT]_y[dOLR]_pole.txt',delimiter=',')
dt1 = data1[:,0]
dolr1 = data1[:,1]

data2 = np.genfromtxt('/Users/nickedkins/Dropbox/Figure Data/x[dT]_y[dOLR]_tropic.txt',delimiter=',')
dt2 = data2[:,0]
dolr2 = data2[:,1]

# plt.plot(dt1,dolr1,'-o',c='b',label='Poles')
# plt.plot(dt2,dolr2,'-o',c='r',label='Tropics')
# plt.xlabel('Change in surface temperature (K)')
# plt.ylabel('Change in OLR (Wm$^{-2}$)')

plt.plot(dolr1,dt1,'-o',c='b',label='Poles')
plt.plot(dolr2,dt2,'-o',c='r',label='Tropics')
plt.ylabel('Change in surface temperature (K)')
plt.xlabel('Change in OLR (Wm$^{-2}$)')

plt.legend()
plt.show()