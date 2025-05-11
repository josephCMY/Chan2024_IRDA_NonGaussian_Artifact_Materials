'''
    SCRIPT TO ILLUSTRATE DRY CLOUD IDEA USING SINGLE OBS ANALYSIS INCREMENTS 

    This script works with obs sequence files from DART.

    Command-line inputs:
    1) Start date in ccyymmddHHMM format
    2) Final date in ccyymmddHHMM format
    3) Time interval



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
date_nw = datetime.strptime('202301130000', '%Y%m%d%H%M')
ens_size=50
targ_latlon = np.array([-10, 82]) 
mem_ind = 0 


'''
    LOAD OBSERVATION PRIORS AND SYNTH OBS DATA
'''
fname = date_nw.strftime('DART_files/NoDA/obs_prior/obs_seq.final.%Y%m%d%H%M')
f = open( fname, 'r' )
obs_dict = read_obs_seq_prior_and_obsval( f, obs_kind=260, ens_size=ens_size )
f.close()



'''
    OBTAINING OBSERVATION AT DESIRED LOCATION
'''
# Determine index corresponding to location
dist = np.sqrt(
    (obs_dict['latitude'] - targ_latlon[0])**2
    + (obs_dict['longitude'] - targ_latlon[1])**2
)
targ_obs_ind = np.argmin(dist)
# print( targ_obs_ind)

# Plot observations
all_obs_prior_names = ['prior ensemble member     %2d' % (ee+1) \
                            for ee in range(ens_size)]
all_obs_priors = np.array(
    [ obs_dict[all_obs_prior_names[ee]][targ_obs_ind] for ee in range(ens_size) ] 
)

obsval = obs_dict['observations'][targ_obs_ind]

print( all_obs_priors)
print( np.mean(all_obs_priors < 230))
print(obsval)
# plt.hist( all_obs_priors, bins=10)
# plt.savefig('tmp.png')

# quit()



'''
    Load NoDA WRF data to compute correlations
'''
# Generate paths to files
NoDA_ens_dir_path = '/fs/ess/PAS2635/Chan2025_IRDA_NonGaussianArtefact_Paper/NoDA/prior/202301130000'
ens_filepath_list = ([
    '/fs/ess/PAS2635/Chan2025_IRDA_NonGaussianArtefact_Paper/NoDA/prior/202301130000/wrfinput_d01_%05d' % (ee+1) 
    for ee in range(ens_size) 
])

# Init variables to hold data
qhydro_vname_list = ['QCLOUD','QICE','QRAIN','QSNOW','QGRAUP']
ens_data_dict = {}
ens_data_dict['qv'] = [None for ee in range(ens_size)]
ens_data_dict['pressure'] = [None for ee in range(ens_size)]
ens_data_dict['t'] = [None for ee in range(ens_size)]
ens_data_dict['qh'] = [None for ee in range(ens_size)]
# ens_data_dict['u'] = [None for ee in range(ens_size)]

# Load data
for ee in range( ens_size ):

    # Open file
    f = ncopen(ens_filepath_list[ee],'r')

    # lat lon arrays
    xlat = f.variables['XLAT'][0,:,0]
    xlon = f.variables['XLONG'][0,0,:]
    base_pres = f.variables['PB'][0]

    # Obtain indices defining zonal vertical cross section in the vicinity of targetted latlon
    lon_ind0 = np.argmin( np.abs( targ_latlon[1] -1 - xlon ) )
    lon_ind1 = np.argmin( np.abs( targ_latlon[1] +1 - xlon ) )
    lat_ind = np.argmin( np.abs( targ_latlon[0] - xlat ) )

    # obtain qvapor, theta and pressure data
    ens_data_dict['qv'][ee] = f.variables['QVAPOR'][0][:,lat_ind, lon_ind0:lon_ind1+1]
    ens_data_dict['t'][ee] = f.variables['T'][0][:,lat_ind, lon_ind0:lon_ind1+1]+300
    ens_data_dict['pressure'][ee] = wrf_getvar(f, 'pres', meta=False, units='Pa')[:,lat_ind, lon_ind0:lon_ind1+1]
    # ens_data_dict['u'][ee] = wrf_getvar(f, 'ua', meta=False)[:,lat_ind, lon_ind0:lon_ind1+1]

    # obtain hydrometeor mixing ratio
    ens_data_dict['qh'][ee] = np.zeros( ens_data_dict['qv'][ee].shape )
    for vname in qhydro_vname_list:
        ens_data_dict['qh'][ee] += f.variables[vname][0][:,lat_ind, lon_ind0:lon_ind1+1]

    # Close file
    f.close()

    print( 'loaded', ee+1 )
# --- end of data reading


# Convert lists to arrays
for vname in ens_data_dict.keys():
    ens_data_dict[vname] = np.array( ens_data_dict[vname] )



'''
    COMPUTE ENSEMBLE CORRELATIONS WITH RESPECT TO BT
'''
ens_corr_dict = {}
irbt_perts = all_obs_priors - np.mean(all_obs_priors)
irbt_sigma = np.std( irbt_perts, ddof=1 )
for vname in ens_data_dict.keys():
    var_perts = ens_data_dict[vname] - np.mean(ens_data_dict[vname], axis=0)
    var_sigma = np.std( ens_data_dict[vname], axis=0, ddof=1 )
    cov = np.zeros_like( var_perts[0] )
    for ee in range( ens_size ):
        cov += (
            (irbt_perts[ee] * var_perts[ee]) 
            / (ens_size-1.)
        )
    ens_corr_dict[vname] = cov / (irbt_sigma * var_sigma)

    # --- end of loop over ensemble
# --- end of loop over variable names



'''
    COMPUTE ANALYSIS INCREMENT AND POSTERIOR FOR MEMBER 2
'''
incre_dict = {}
poste_dict = {}
obs_err_sigma = 1.
for vname in ens_data_dict.keys():
    # Generate covariance
    var_perts = ens_data_dict[vname] - np.mean(ens_data_dict[vname], axis=0)
    cov = np.zeros_like( var_perts[0] )
    for ee in range( ens_size ):
        cov += (
            (irbt_perts[ee] * var_perts[ee]) 
            / (ens_size-1.)
        )

    # Kalman gain and square root factor
    kalman_gain = cov / (irbt_sigma**2 + obs_err_sigma**2)
    phi = 1./( 
        1 + np.sqrt( 
            obs_err_sigma**2
            / (irbt_sigma**2 + obs_err_sigma**2)
        )
    )

    # innovation
    inno = obsval - np.mean( all_obs_priors )
    # print( inno )

    # obs perturbation
    yprime = all_obs_priors[mem_ind] - np.mean( all_obs_priors )

    # increment calculation
    incre_dict[vname] = kalman_gain * ( inno - phi*yprime)
    poste_dict[vname] = ens_data_dict[vname][mem_ind]  + incre_dict[vname]

# --- end of loop over variable names



'''
    Function to compute saturation ratios
'''
def calc_saturation_ratio(qvapor, theta, pressure ):

    # Convert theta to temperature in kelvins
    tk = theta * np.power( pressure / 1e5, 0.286 )

    # compute vapor pressure
    e = ( qvapor/ (287/461 +qvapor) ) * pressure

    # compute saturation vapor pressures
    e_sat = np.empty_like(e)
    flag_liq = (tk >= 273.16)
    flag_ice = (tk <  273.16)
    e_sat[flag_liq] = (
        611 
        * np.exp( 
            6808 * (1/273.16 - 1/tk[flag_liq])
            - 5.09 * np.log( tk[flag_liq]/273.16 )
        )
    )
    e_sat[flag_ice] = (
        611 
        * np.exp( 
            6293 * (1/273.16 - 1/tk[flag_ice])
            - 0.555 * np.log( tk[flag_ice]/273.16 )
        )
    )

    # Return saturation ratios
    return e/e_sat



'''
    GENERATE DICTIONARY TO HOLD PRIOR MEMBER mem_ind
'''
prior_dict = {}
for vname in ens_data_dict.keys():
    prior_dict[vname] = ens_data_dict[vname][mem_ind]




'''
    COMPUTE SATURATION RATIOS BEFORE AND AFTER DA
    Note: Prior pressure used in calculation because not updated by DART.
'''
for dicttype in [prior_dict, poste_dict]:
    dicttype['saturation'] = calc_saturation_ratio( dicttype['qv'], dicttype['t'], prior_dict['pressure'] )
    # print( dicttype['saturation'].min(), dicttype['saturation'].max())
    # dicttype['pressure'] -= base_pres



'''
    CONVENIENT CHANGE OF UNITS -- kg/kg to g/kg
'''
for dicttype in [prior_dict, poste_dict, incre_dict]:
    for vname in ['qh','qv']:
        dicttype[vname] *= 1e3




'''
    VISUALIZE!
    Row 1: Forecast and analysis qhydro and S
    Row 2: Correlations
    Row 3: Analysis increments
'''
fig, axs = plt.subplots( nrows=3, ncols=2, figsize = (7,6))

# Color ranges
crange_dict = {
    'qv': np.linspace(0,18, 7), #np.power(10, np.linspace(-4,2, 11)),
    't': np.linspace(290,390,6),
    'qh': np.linspace(0.01,1.01, 6),
    'u': np.linspace(-10,10,9),
}

incre_crange_dict = {
    'qv': np.linspace(-1.2,1.2, 7), #np.power(10, np.linspace(-4,2, 11)),
    'qh': np.linspace(-1.2,1.2, 7),
}

units_dict = {
    'qv': 'g/kg',
    'qh': 'g/kg'
}

cmap_dict = {
    'qv': 'Greens',
    'qh': 'inferno_r',
    'u': 'PRGn',
    't': 'plasma'
}

# lon mesh
lonmesh, tmp = np.meshgrid( xlon[lon_ind0:lon_ind1+1], np.arange(44) )



'''
    Row 1: Prior and posterior clouds
'''
vname = 'qh'
# Prior clouds
cnf = axs[0,0].contourf( 
        lonmesh, prior_dict['pressure']/100., prior_dict[vname]+1e-7, 
        crange_dict[vname], cmap = cmap_dict[vname],
        extend='max'
)
plt.colorbar(cnf, ax=axs[0,0], label='%s (%s)' % (vname, units_dict[vname]))
axs[0,0].set_title('a%d) Fcst Mem %d %s ' % (1, mem_ind+1, vname.upper()), loc='left' )


# Poste clouds
cnf = axs[0,1].contourf( 
        lonmesh, poste_dict['pressure']/100., poste_dict[vname]+1e-7, 
        crange_dict[vname], cmap = cmap_dict[vname],
        extend='max'
)
plt.colorbar(cnf, ax=axs[0,1], label='%s (%s)' % (vname, units_dict[vname]))
axs[0,1].set_title('a%d) Anlys Mem %d %s ' % (2, mem_ind+1, vname.upper()), loc='left' )




'''
    Other rows
'''
# Loop over all variable names
for icol, vname in enumerate( ['qv','qh'] ):

    # Plot correlations
    cnf = axs[1, icol].contourf( 
        lonmesh, prior_dict['pressure']/100., ens_corr_dict[vname], 
        np.linspace(-0.9,0.9,7), cmap = 'RdBu_r'
    )
    plt.colorbar(cnf, ax=axs[1, icol], orientation='vertical', label='Corr (no units)')
    axs[1, icol].set_title('b%d) Corr(IR BT, %s)' % (icol+1, vname.upper()), loc='left')

    # Plot increment
    absvals = np.abs( incre_dict[vname])
    absvals = absvals[ (np.isinf(absvals) + np.isnan(absvals) == 0) ]
    limit = np.percentile( absvals, 95)
    crange = np.linspace( limit*-1, limit, 11 )
    cnf = axs[2, icol].contourf( 
        lonmesh, prior_dict['pressure']/100., incre_dict[vname], 
        incre_crange_dict[vname], cmap = 'PuOr_r', extend='both'
    )
    plt.colorbar(cnf, ax=axs[2, icol], orientation='vertical', label='%s incre (%s)' % (vname, units_dict[vname]))
    axs[2, icol].set_title('c%d) Mem %d Incre %s' % (icol+1, mem_ind+1, vname.upper()), loc='left')


# --- end of loop over rows



# overlay hatches in prior plot to indicate saturation
cn1 = axs[0,0].contourf( 
    lonmesh, prior_dict['pressure']/100., prior_dict['saturation'], 
    [0.9, 10], hatches = ['xx'], colors='none'
)
axs[0,0].contour( 
    lonmesh, prior_dict['pressure']/100., prior_dict['saturation'], 
    [0.9], colors = 'dodgerblue', linewidths=1, linestyles = ['-']
)
for collection in cn1.collections:
    collection.set_edgecolor( 'dodgerblue' )
    

# overlay hatches in posterior plot to indicate saturation
cn1 = axs[0,1].contourf( 
    lonmesh, prior_dict['pressure']/100., poste_dict['saturation'], 
    [0.9, 10], hatches = ['xx'], colors='none'
)
axs[0,1].contour( 
    lonmesh, prior_dict['pressure']/100., poste_dict['saturation'], 
    [0.9], colors = 'dodgerblue', linewidths=1, linestyles = ['-']
)
for collection in cn1.collections:
    collection.set_edgecolor( 'dodgerblue' )


for irow in range(3):
    for icol in range(2):
        axs[irow, icol].set_ylim([999,99])
        axs[irow, icol].set_ylabel('Pres (hPa)')
        axs[irow, icol].set_xlabel('Lon (deg E)')
        axs[irow, icol].axvline(targ_latlon[1], linestyle='--', color='k')

plt.suptitle( 
    'Single Obs Test (Obs: %d K, Mem %d: %d K)' 
    % (obsval, mem_ind+1, all_obs_priors[mem_ind]) 
)

plt.tight_layout()
plt.savefig('figs/single_obs_test.pdf')