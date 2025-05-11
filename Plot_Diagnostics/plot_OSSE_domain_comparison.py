'''
    SCRIPT TO COMPARE SOME DOMAIN LEVELS ACROSS EXPERIMENTS


'''

from netCDF4 import Dataset as ncopen
import numpy as np
from matplotlib import use
use('agg')
import matplotlib.pyplot as plt
# import cartopy.crs as ccrs
# from sys import argv
from datetime import datetime, timedelta
from copy import deepcopy
import matplotlib.colors as mcolors
from math import pi as PI
import pickle
from gc import collect as gc_collect

import matplotlib.dates as mdates
myFmt = mdates.DateFormatter('%b-%d\n%HUTC')
myBlankFmt = mdates.DateFormatter(' ')

from pylib_expt_specific_settings import *
from pylib_read_obs_priors import *

from wrf import getvar as wrf_getvar

# from numba import njit

# import matplotlib.patches as patches





'''
    USER INPUTS
'''
date_nw = datetime.strptime('202301151200', '%Y%m%d%H%M')
ens_size=50
expt_list = ['NoDA','BASIC','FCSTCLOUD']
vname_list = ['QVAPOR','T','ua','va']
flag_compute_from_raw_data = False
pkl_fname = 'pkl_files/domain_comparison.pkl'

# Plotting model levels
plot_lvls = [10, 24]


'''
    USING HYPOTHETICAL PRESSURE LEVELS
'''

eta_lvls =  np.array([ 0.996521  , 0.9881265 , 0.976891 ,  0.962891  , 0.946225  , 0.927016,
                       0.9054105 , 0.88156605, 0.8556515,  0.8278555 , 0.79838955, 0.7674775,
                       0.73533046, 0.70214903, 0.6681355,  0.63348997, 0.59841347, 0.5631155,
                       0.5278065 , 0.49267948, 0.457909 ,  0.4236615 , 0.390096  , 0.357357,
                       0.3255645 , 0.2948205 , 0.265225 ,  0.2368785 , 0.2098805 , 0.1843245,
                       0.160296  , 0.1378805 , 0.1171655,  0.09825251, 0.0812925 , 0.0664215,
                       0.0535585 , 0.0424325 , 0.0328125,  0.0245355 , 0.017449  , 0.0114115,
                       0.0062925 , 0.0019735 ] )
plvls = eta_lvls* (1e5-2e3) +2e3
plvls /= 1e2    



'''
    Load ensemble WRF data and nature run. 
'''
# Will only run this calculation if needed. Very time consuming
if flag_compute_from_raw_data:

    # Prep dictionary to hold data
    data_dict = {}
    data_dict['Nature Run'] = {}
    for ename in expt_list:
        data_dict[ename] = {}
        for vname in vname_list:
            data_dict[ename][vname] = None


    # Loop over experiments
    for ename in expt_list:
        
        # Ens dir path
        ens_dir_path = ( 
            '/fs/ess/PAS2635/Chan2025_IRDA_NonGaussianArtefact_Paper/%s/prior/%s'
            % (ename, date_nw.strftime('%Y%m%d%H%M'))
        )

        # Load data from each member
        for ee in range(ens_size):

            print('Loading %s member %04d' % (ename, ee+1))

            # File opening
            fname = '%s/wrfinput_d01_%05d' % (ens_dir_path, ee+1)
            f = ncopen( fname, 'r')

            # Loading and averaging variable across ens
            for vname in vname_list:

                # Load and average
                v_dom = wrf_getvar( f, vname, meta=False )
                v_dom /= ens_size

                # Special treatment for first member
                if ee == 0:
                    data_dict[ename][vname] = np.zeros_like( v_dom )
                
                # Increment data into dictionary
                data_dict[ename][vname] += v_dom
            
            # --- End of loop over variables

            # Close file to release handle
            f.close()

            # Garbage collection to speed up code
            if ee % 10 == 0 and ee > 0:
                gc_collect()

        # --- End of loop over ensemble

        # Garbage collection to speed up code
        gc_collect()
    # --- End of loop over OSSEs


    # Load from nature run
    # --------------------
    ename = 'Nature Run'

    # file reading
    fname = date_nw.strftime( 'nature_run/wrfinput_d01_%Y-%m-%d_%H:%M:%S' )
    f  = ncopen( fname, 'r')

    # Load and average
    for vname in vname_list:
        data_dict[ename][vname] = wrf_getvar( f, vname, meta=False )
    # -- end of loop over variables

    data_dict['lat'] = np.squeeze( f.variables['XLAT'])
    data_dict['lon'] = np.squeeze( f.variables['XLONG'])


    # Store into pickle file
    with open( pkl_fname, 'wb') as f:
        pickle.dump( data_dict, f )
    # -- end of pickling

# Handling situation where the calcualtion flag is false.
else:
    with open( pkl_fname, 'rb') as f:
            data_dict= pickle.load( f )




'''
    PLOTTING TIME!
'''
fig, axs = plt.subplots( nrows=4, ncols=2, figsize=(6,8))

lonmesh, latmesh = data_dict['lon'][::40,::40], data_dict['lat'][::40,::40]
# Generate plots of domain
for ip, p_ind in enumerate( plot_lvls ): 

    t_crange = [
        np.linspace(0,5, 6)+4 ,
        np.linspace(0,5, 6)+37
    ][ip]

    for ie, ename in enumerate( ['Nature Run', 'NoDA','BASIC','FCSTCLOUD'] ):

        # Select axis
        ax = axs[ie, ip]

        # # Plot temperature field
        # cnf = ax.contourf( 
        #     data_dict['lon'], data_dict['lat'],
        #     data_dict[ename]['T'][p_ind], t_crange,
        #     cmap = 'RdBu_r', extend='both'
        # )
        # plt.colorbar( cnf, ax=ax)

        # Plot flows
        sp = ax.quiver(
            lonmesh, latmesh,
            data_dict[ename]['ua'][p_ind,::40,::40], data_dict[ename]['va'][p_ind,::40,::40],
        )
        
        # # Plot flow patterns
        # sp = ax.streamplot(
        #     data_dict['lat'][::10], plvls[::-1],
        #     data_dict[ename]['V'][::-1,::10], data_dict[ename]['omega'][::-1,::10]*-1,
        #     broken_streamlines=False
        # )

        # Labelling
        ax.set_title(ename)
        # ax.set_ylim([1000,0])

plt.savefig('tmp.png')

    
