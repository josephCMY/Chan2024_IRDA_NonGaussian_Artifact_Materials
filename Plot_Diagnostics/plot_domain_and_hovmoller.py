'''
    SCRIPT TO COMPARE NATURE RUN WV-BT AGAINST OBSERVED WV-BT

    This script works with obs sequence files from DART.

    Command-line inputs:
    1) Start date in ccyymmddHHMM format
    2) Final date in ccyymmddHHMM format
    3) Time interval

    Will compare WV-BT channel observations

    Assumes:
    1) The observation values live in osse_met9_seviri_channel08_obs_seq/ccyymmddHHMM/obs_seq.in
    2) The simulated obs values live in osse_met9_seviri_channel08_obs_seq/ccyymmddHHMM/obs_seq.out
    3) Both obs_seq.in and obs_seq.out only contain SEVIRI channel 8 SEVIRI values

'''

from netCDF4 import Dataset as ncopen
import numpy as np
from matplotlib import use
use('agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from sys import argv
from datetime import datetime, timedelta
from copy import deepcopy
import matplotlib.colors as mcolors
from math import pi as PI

import matplotlib.dates as mdates
myFmt = mdates.DateFormatter('%b-%d\n%HUTC')
myBlankFmt = mdates.DateFormatter(' ')

from pylib_expt_specific_settings import *

from numba import njit

import matplotlib.patches as patches








''' Function to compute 11.2 um Window-bT from OLR'''
# (1984) as given in Yang and Slingo (2001)
# Tf = tb(a+b*Tb) where a = 1.228 and b = -1.106e-3 K^-1
# OLR = sigma*Tf^4 
# where sigma = Stefan-Boltzmann constant = 5.67x10^-8 W m^-2 K^-4
def olr2bt( OLR ):
    a = 1.228
    b = -1.106e-3
    sigma = 5.67e-8 # W m^-2 K^-4
    tf = (OLR/sigma)**0.25
    tb = (-a + np.sqrt(a**2 + 4*b*tf))/(2*b)
    return tb









'''
    FUNCTION TO EXTRACT AN OBSERVATION TYPE FROM DART OBSERVATION FILE OBS_SEQ.FINAL
'''
def read_obs_seq( obs_seq_file_handle, obs_kind=260, entryname='observation' ):
    
    # Load ASCII data
    all_lines = obs_seq_file_handle.readlines()


    # Determine number of desired observations by reading header info
    # ----------------------------------------------------------------
    
    # OBS_KIND_SEARCH LOOP: Looking for header line with correct obs kind code
    for ll in range(len(all_lines)):

        # Load string corresponding to line ll
        line0 = all_lines[ll]

        # Check the first 12 characters for obs_kind code
        obs_kind_str = '%d' % obs_kind
        if ( obs_kind_str in line0[:12] ):
            # If the correct obs_kind code is found, break the search loop
            break
        # --- End of if statement

    # --- End of OBS_KIND_SEARCH LOOP

    # Read number of observations 
    line_with_obs_count = all_lines[ll+2]
    num_obs = int( deepcopy(line_with_obs_count).split()[1] )
        
    # ---------------- End of determining number of desired observations





    # Determine ordering of observation entry (obs, truth, prior mean, etc)
    # in every OBS entry
    # Want: observations, prior ensemble mean, prior ensemble spread,
    #       posterior ensemble mean, and posterior ensemble spread
    # Note that the indices are relative to the line that says "OBS"
    # ----------------------------------------------------------------------------
    ind_obs_entry_dict = {}
    entries_wanted = [entryname]

    # Detect where 'observations' key occurs
    for ll in range( len(all_lines) ):
        if ( 'observation' in all_lines[ll] ):
            ref_ll = deepcopy( ll )
            break
        # --- End of if statement to find string that says 'observations'
    # --- End of loop search for the string 'observations'

    # PARSE_OBS_ENTRY_ORDERING_HEADER LOOP: Loop to determine how different obs values are
    # ordered
    for ll in range( ref_ll, len(all_lines) ):
        # Loop over all entries desired
        for entry in entries_wanted:
            if ( entry in all_lines[ll] ): 
                ind_obs_entry_dict[entry] = ll - ref_ll + 1
            # --- End of checking for entry
        # End of loop over all entries wanted
    
        # If all desired entries have been detected, stop PARSE_OBS_ENTRY_ORDERING_HEADER
        if ( len( ind_obs_entry_dict.keys() ) ==  len(entries_wanted) ):
            break
        # --- End of checking whether all desired entries are detected

    # --- End of PARSE_OBS_ENTRY_ORDERING_HEADER LOOP

    # ------------------------- Finished parsing ordering of observation values

    

    # Init dictionary to hold all observation entries desired
    # --------------------------------------------------------
    obs_dict = {}
    obs_dict['longitude']               = np.zeros( num_obs, dtype=float ) +np.nan
    obs_dict['latitude']                = np.zeros( num_obs, dtype=float ) +np.nan
    obs_dict['level']                   = np.zeros( num_obs, dtype=float ) +np.nan
    obs_dict['vertical coord type']     = np.zeros( num_obs, dtype=int ) +np.nan
    for entry in entries_wanted:
        obs_dict[entry] = np.zeros( num_obs, dtype=float ) +np.nan
    # --------------------------------- End of init obs entry holder
    
    
    
    # Detect lines where 'OBS' string occurs
    obs_chunk_start_line_numbers = np.zeros(999999, dtype=int) - 999
    obs_cnt = 0
    for ll in range( len(all_lines) ):

        # If "OBS" occurs, then mark this line index
        if ("OBS  " in all_lines[ll] ):
            obs_chunk_start_line_numbers[obs_cnt] = ll
            obs_cnt += 1
        # --- End of if statment

        # Special safety statement
        if obs_cnt >= len(obs_chunk_start_line_numbers):
            print('Ran out of line numbers!')
            print('Please increase size of obs_chunk_start_line_numbers')
            print('Exiting')
            quit()
        # --- End of safety statement
    
    # --- End of loop to search for "OBS"



    # Look for desired observation kind and load data
    # -----------------------------------------------

    desired_obs_cnt = 0

    # SEARCH_ALL_OBS_CHUNKS LOOP: Search over all obs chunks
    for obs_ind in range(obs_cnt):

        # Determine range of lines containing the observation chunk
        ll0 = obs_chunk_start_line_numbers[obs_ind]
        if ( obs_chunk_start_line_numbers[obs_ind+1] < 0  ):
            obs_string_chunk = all_lines[ll0:]
        else:
            ll1 = obs_chunk_start_line_numbers[obs_ind+1]
            obs_string_chunk = all_lines[ll0:ll1]
        # --- End of observation chunk determination


        # Determine this OBS's obs_kind code
        for ll in range(len(obs_string_chunk)):
            # Checking if 'obs_kind' string detected
            if ( 'kind' in obs_string_chunk[ll] ):
                current_obs_kind = int(obs_string_chunk[ll+1])
            # --- End of checking for 'obs_kind' string
        # --- End of loop over lines in obs_string_chunk


        # Skip this OBS chunk the obs_kind for this OBS is not the desired one
        if ( current_obs_kind != obs_kind ):
            continue
        # --- End of skipping statement


        # Read in observation entries
        for entry in entries_wanted:
            relative_index = ind_obs_entry_dict[entry]
            obs_dict[entry][desired_obs_cnt] = float( obs_string_chunk[relative_index])
        # --- End of obs entry reading


        # SEARCH_LOCATION LOOP: Search thru obs_string_chunk for obs location info
        for ll in range(len(obs_string_chunk)):

            # Check if 'loc3d' header is detected
            if ( 'loc3d' in obs_string_chunk[ll] ):
                # Read the next line for location values
                location_vals = (obs_string_chunk[ll+1]).split()
                obs_dict['longitude'][desired_obs_cnt]           = float( location_vals[0] ) * 180/ PI
                obs_dict['latitude'][desired_obs_cnt]            = float( location_vals[1] ) * 180/ PI
                obs_dict['level'][desired_obs_cnt]               = float( location_vals[2] )
                obs_dict['vertical coord type'][desired_obs_cnt] =   int( location_vals[3] )

                # Exit SEARCH_LOCATION loop
                break
            
            # --- End of if-statement for detecting 'loc3d' header

        # --- End of SEARCH_LOCATION LOOP

        # Since we have found a desired obs, increment desired_obs_cnt 
        desired_obs_cnt += 1

    # --- End of SEARCH_ALL_OBS_CHUNKS LOOP

    # ---------------------------------------- Finished searching through 
    #                                          all obs chunks for desired obs


    # Check if we have the correct number of desired observations loaded
    if ( desired_obs_cnt != num_obs ):
        print('OH NO! Number of loaded observations do not match the number of observations detected')
        print('Exiting!')
        quit()
    # --- End of safety check

    # Return all desired observations
    return obs_dict











'''
    FUNCTION TO GENERATE LON-TIME (HOVMOLLER) DIAGRAM ARRAY FROM OBSERVATIONS

    Inputs:
    1) obs_dict_list -- List of obs dicts (1 per date)
'''
def bin_obs_into_hovmoller_array( obs_dict_list, entryname = "observation", lon_bin_edges = np.linspace(50,90,41)*1. ):

    # Determine number of dates & lon bins
    nDates = len(obs_dict_list)
    nLon = len(lon_bin_edges[1:])

    # Init array for hovmoller diagram
    hov_array = np.zeros( [nDates, nLon] )
    
    # Generate hovmoller array
    for dd in range(nDates):

        # Isolate only 5N to 15S
        mask = ( ( (obs_dict_list[dd]['latitude'] < 5) * (obs_dict_list[dd]['latitude'] > -15) ) > 0.5 )
        obsvals = obs_dict_list[dd][entryname][mask]
        obslons = obs_dict_list[dd]['longitude'][mask]

        hov_array[dd,:] = compute_avg_in_lon_bins( 
            deepcopy( obsvals ),
            deepcopy( obslons ),
            lon_bin_edges
        )

    # Output bin edge array and hovmoller array

    return lon_bin_edges, hov_array

    


# Quick function to compute averages over lon bins
@njit
def compute_avg_in_lon_bins(obs_vals, obs_lon, lon_bin_edges):

    # Determine number of bins
    nBins = len(lon_bin_edges[1:])

    # Determine number of obs
    nObs = len(obs_vals)

    # Array to track number of obs in each bin
    cnt_arr = lon_bin_edges[1:]*0

    # init array to hold avg values
    avg_arr = lon_bin_edges[1:]*0

    # Loop over all obs
    for iob in range(nObs):

        # Loop over all bins
        for bb in range(nBins):

            # Check if obs is within bin
            if ( lon_bin_edges[bb] <= obs_lon[iob] 
                 and obs_lon[iob] < lon_bin_edges[bb+1] ):
                
                # If obs is within bin, include and escape loop
                cnt_arr[bb] += 1
                avg_arr[bb] += obs_vals[iob]
                break
            
            # --- End of if statement
        # --- End of loop over all bins
    # --- End of loop over all obs

    # Compute average and return
    avg_arr /= cnt_arr

    return avg_arr

        






















'''
    LOAD USER INPUTS
'''
date_st = datetime.strptime(argv[1], '%Y%m%d%H%M')
date_ed = datetime.strptime(argv[2], '%Y%m%d%H%M')
time_interval = float( argv[3] )











'''
    DETERMINE LIST OF DATES IN THE DATE RANGE
'''
# Construct list of dates
date_list = []
date_nw = deepcopy( date_st ) 
while ( date_nw <= date_ed ):
    date_list.append( date_nw )
    date_nw += timedelta( minutes = time_interval )

# Determine number of dates
nDates = len( date_list )


dateticks = [ date_list[0], date_list[0]+ (date_list[-1]-date_list[0])/3, 
              date_list[0]+ (date_list[-1]-date_list[0])*2/3, \
              date_list[-1] ]





'''
    LOAD REAL OBSERVATION DATA
'''
# Init array to hold all data entries
all_real_obs = [None for i in range(nDates)]

# Loop over all dates
for dd in range(nDates):

    date_nw = date_list[dd]

    # Load real data 
    fname = date_nw.strftime('OSSE_obs_seq/%Y%m%d%H%M/obs_seq.in')
    f = open( fname, 'r' )
    all_real_obs[dd] = read_obs_seq( f, obs_kind=260, entryname='observation' )
    f.close()

# --- End of loop over all dates












'''
    LOAD NATURE RUN OBSERVATION DATA
'''
# Init array to hold all data entries
all_nature_obs = [None for i in range(nDates)]

# Loop over all dates
for dd in range(nDates):

    date_nw = date_list[dd]

    # Load real data 
    fname = date_nw.strftime('OSSE_obs_seq/%Y%m%d%H%M/obs_seq.out')
    f = open( fname, 'r' )
    all_nature_obs[dd] = read_obs_seq( f, obs_kind=260, entryname='truth' )
    f.close()

# --- End of loop over all dates




# '''
#     Load OLR
# '''
# # Load OLR
# f = ncopen( 'nature_run_olr.nc', 'r')
# olr = f.variables['OLR'][5]
# times = f.variables['Times'][5]
# f.close()

# windowbt = olr2bt( olr )





'''
    GENERATE HOVMOLLER ARRAY FOR REAL OBSERVATION DATA & NATURE RUN OBSERVATION DATA
'''

lon_bin_edges, nature_hov_arr = bin_obs_into_hovmoller_array( all_nature_obs, entryname = "truth" )
lon_bin_edges, real_hov_arr   = bin_obs_into_hovmoller_array( all_real_obs  , entryname = "observation" )



'''
    QUANTILE-MATCH THE HOVMOLLER ARRAYS TO SIMPLIFY VISUALIZATION
'''

# Determine ranks
shp = nature_hov_arr.shape
rank_list = np.arange(np.prod(shp))
nature_hov_rank = (nature_hov_arr.argsort(axis=None)).argsort(axis=None).reshape(shp) *100./np.prod(shp)
real_hov_rank = (real_hov_arr.argsort(axis=None)).argsort(axis=None).reshape(shp) *100./np.prod(shp)








'''
    LOAD COORDINATES OF DOMAIN
'''
# Load coordinate
f = ncopen('nature_wrfinput_d01_202301130000.nc','r')
xlat = f.variables['XLAT'][0]
xlon = f.variables['XLONG'][0]












'''
    PLOT OUT DOMAIN AND HOVMOLLER DIAGRAMS
'''

# Init figure
fig = plt.figure( figsize=(9,6))
axs = [None, None, None, None]

# Figure axis for domain plot
axs[0] = plt.axes([0.05, 0.61, 0.4, 0.35 ], projection=ccrs.PlateCarree())

# Figure axis for domain plot
axs[1] = plt.axes([0.05, 0.18, 0.4, 0.35 ], projection=ccrs.PlateCarree())

# Figure axis for hovmoller plots
axs[2] = fig.add_axes( [ 0.55, 0.2, 0.2, 0.75 ] )
axs[3] = fig.add_axes( [ 0.78, 0.2, 0.2, 0.75 ] )



# # Init figure
# fig = plt.figure( figsize=(6,9))
# axs = [None, None, None, None]

# # Figure axis for domain plot
# axs[0] = plt.axes([0.1, 0.63, 0.4, 0.32 ], projection=ccrs.PlateCarree())

# # Figure axis for domain plot
# axs[1] = plt.axes([0.55, 0.63, 0.4, 0.32 ], projection=ccrs.PlateCarree())

# # Figure axis for hovmoller plots
# axs[2] = fig.add_axes( [0.1, 0.15, 0.4, 0.3 ] )
# axs[3] = fig.add_axes( [ 0.55, 0.15, 0.4, 0.3 ] )






# Plot nature run domain
ax = axs[0]
img = plt.imread('bluemarble-2048.png')
img_extent = (-180, 180, -90, 90)
ax.imshow(img, origin='upper', extent=img_extent, transform=ccrs.PlateCarree() )
# ax.coastlines(resolution='50m')

# Mark domain boundaries
rect = patches.Rectangle(
        (xlon.min(), xlat.min()), 
        xlon.max()-xlon.min(), 
        xlat.max()-xlat.min(), linewidth=1, 
        edgecolor='none', facecolor="white", alpha = 0.3) 
ax.add_patch(rect) 
ax.plot( xlon[0,:], xlat[0,:], linewidth=6, color='r', zorder=10)
ax.plot( xlon[-1,:], xlat[-1,:], linewidth=6, color='r', zorder=10)
ax.plot( xlon[:,0], xlat[:,0], linewidth=6, color='r', zorder=10)
ax.plot( xlon[:,-1:], xlat[:,-1], linewidth=6, color='r', zorder=10)


gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
gl.top_labels = False
gl.right_labels = False
ax.set_extent([40, 100, -26, 6], crs=ccrs.PlateCarree())


dd=0
cnf = ax.tricontourf( all_nature_obs[dd]['longitude'], all_nature_obs[dd]['latitude'],
                  all_nature_obs[dd]['truth'], np.linspace(204,224,6),
                  cmap = 'inferno_r', extend='min', zorder=3
             )

ax.tricontourf( all_nature_obs[dd]['longitude'], all_nature_obs[dd]['latitude'],
                  all_nature_obs[dd]['truth'], np.linspace(200,234,2),
                  cmap = 'Blues', extend='min', zorder=2
             )

# ax.contourf( xlon, xlat, windowbt, 
#                  [100,280], cmap = 'Blues', zorder=2
#              )

ax.set_title( date_list[dd].strftime( 'a) Nature Run at %H:%M UTC on %d-%b-%Y'), loc='left')

# cbar = fig.colorbar(cnf, ax=ax, orientation='vertical')
# cbar.ax.set_ylabel('SEVIRI 6.2 um BT (K)')









# Plot real observation over domain
ax = axs[1]
img = plt.imread('bluemarble-2048.png')
img_extent = (-180, 180, -90, 90)
ax.imshow(img, origin='upper', extent=img_extent, transform=ccrs.PlateCarree() )
# ax.coastlines(resolution='50m')

# Mark domain boundaries
rect = patches.Rectangle(
        (xlon.min(), xlat.min()), 
        xlon.max()-xlon.min(), 
        xlat.max()-xlat.min(), linewidth=1, 
        edgecolor='none', facecolor="white", alpha = 0.3) 
ax.add_patch(rect) 
ax.plot( xlon[0,:], xlat[0,:], linewidth=6, color='r', zorder=10)
ax.plot( xlon[-1,:], xlat[-1,:], linewidth=6, color='r', zorder=10)
ax.plot( xlon[:,0], xlat[:,0], linewidth=6, color='r', zorder=10)
ax.plot( xlon[:,-1:], xlat[:,-1], linewidth=6, color='r', zorder=10)


gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
gl.top_labels = False
gl.right_labels = False
ax.set_extent([40, 100, -26, 6], crs=ccrs.PlateCarree())


dd=0
cnf_domain = ax.tricontourf( all_real_obs[dd]['longitude'], all_real_obs[dd]['latitude'],
                  all_real_obs[dd]['observation'], np.linspace(204,224,6),
                  cmap = 'inferno_r', extend='min', zorder=3
             )

ax.tricontourf( all_real_obs[dd]['longitude'], all_real_obs[dd]['latitude'],
                  all_real_obs[dd]['observation'], np.linspace(200,234,2),
                  cmap = 'Blues', extend='min', zorder=2
             )

# ax.contourf( xlon, xlat, windowbt, 
#                  [100,280], cmap = 'Blues', zorder=2
#              )

ax.set_title( date_list[dd].strftime( 'b) Real Obs at %H:%M UTC on %d-%b-%Y'), loc='left')

# cbar = fig.colorbar(cnf, ax=ax, orientation='vertical')
# cbar.ax.set_ylabel('SEVIRI 6.2 um BT (K)')














# Plot osse data
lon_bin_ctrs = (lon_bin_edges[1:]+lon_bin_edges[:-1])/2
cnf= axs[2].contourf(lon_bin_ctrs, date_list, nature_hov_rank, np.linspace(0,100,11), cmap='RdBu_r')
axs[2].set_title('c) Nature Run', loc='left')



# Plot real data
cnf= axs[3].contourf(lon_bin_ctrs, date_list, real_hov_rank, np.linspace(0,100,11), cmap='RdBu_r')
axs[3].set_title('d) Real Obs', loc='left')

print( np.corrcoef( real_hov_rank.flatten(), y=nature_hov_rank.flatten()  ) )


# Some cosmetic adjustments
axs[2].yaxis.set_major_formatter(myFmt)
axs[3].yaxis.set_major_formatter(myBlankFmt)
axs[2].set_xlabel('Longitude (degrees E)')
axs[3].set_xlabel('Longitude (degrees E)')
axs[2].set_yticks(dateticks)
axs[3].set_yticks(dateticks)

# fig.subplots_adjust(left=0.1, right=0.98, top=0.93, bottom=0.3, wspace=0.2)

# Colorbar for the IR BT domain plots
cbar_ax = fig.add_axes([0.1, 0.075, 0.3, 0.03])
cbar = fig.colorbar(cnf_domain, cax=cbar_ax, orientation='horizontal')
cbar.ax.set_xlabel('SEVIRI 6.2 um BT (K)')

# Colorbar for IR BT percentile plots
cbar_ax = fig.add_axes([0.62, 0.075, 0.3, 0.03])
fig.colorbar(cnf, cax=cbar_ax, orientation='horizontal')
cbar_ax.set_xlabel('SEVIRI 6.25 um Hovmoller Percentiles')





plt.savefig('figs/nature_run_hovmoller.svg')