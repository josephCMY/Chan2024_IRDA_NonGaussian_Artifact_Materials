'''
    SCRIPT TO COMPARE ENSEMBLE MEAN SEVIRI CHANNEL 08 VALUES AGAINST THE OBSERVED VALUES

    This script works with obs_seq.final obs sequence files from DART.

    Command-line inputs:
    1) Start date in ccyymmddHHMM format
    2) Final date in ccyymmddHHMM format
    3) Time interval between plots in minutes

    Will compare OLR-based channel 13, cloud-top temperature and potential vorticity 
    from WRF files against window channel observations (SEVIRI channel 13)
    This is just for qualitative comparison purposes.

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

from pylib_expt_specific_settings import *















'''
    FUNCTION TO EXTRACT AN OBSERVATION TYPE FROM DART OBSERVATION FILE OBS_SEQ.FINAL
'''
def read_obs_seq_final_obs_and_prior( obs_seq_file_handle, obs_kind=260 ):
    
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
    entries_wanted = ['observations', 'prior ensemble mean', 'prior ensemble spread', \
                       'posterior ensemble mean', 'posterior ensemble spread']

    # Detect where 'observations' key occurs
    for ll in range( len(all_lines) ):
        if ( 'observations' in all_lines[ll] ):
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
    LOAD USER INPUTS
'''
date_st = datetime.strptime(argv[1], '%Y%m%d%H%M')
date_ed = datetime.strptime(argv[2], '%Y%m%d%H%M')
time_interval = float( argv[3] )
expt_name = argv[4]
obs_err_sigma=float(argv[5])











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











'''
    LOAD OSSE TRUE OBSERVATION DATA
'''
# Init array to hold all data entries
all_obs_seq_final = [None for i in range(nDates)]


# Loop over all dates
for dd in range(nDates):

    date_nw = date_list[dd]

    # Load real data 
    fname = date_nw.strftime('DART_files/'+expt_name+'/obs_seq.final.%Y%m%d%H%M')
    f = open( fname, 'r' )
    all_obs_seq_final[dd] = read_obs_seq_final_obs_and_prior( f, obs_kind=260 )
    f.close()

# --- End of loop over all dates






'''
    PLOT OUT TIME-SERIES OF BIASES, RMSES AND SPREAD
'''


# Init arrays to hold metrics
omb_bias = np.zeros(nDates)
omb_rmse  = np.zeros(nDates)
omb_rmss = np.zeros(nDates)
oma_rmse = np.zeros(nDates)
oma_bias = np.zeros(nDates)

# Compute metrics for all dates
for dd in range(nDates):

    # Removing strange values...
    flags = (all_obs_seq_final[dd]['posterior ensemble mean'] < 100) + (all_obs_seq_final[dd]['posterior ensemble mean'] > 350 )
    flags += (all_obs_seq_final[dd]['observations'] < 100) + (all_obs_seq_final[dd]['observations'] > 350 )
    flags = (flags == 0)
    # print( np.sum(flags), len(flags))

    omb = (all_obs_seq_final[dd]['observations'] - all_obs_seq_final[dd]['prior ensemble mean'])[flags]
    omb_bias[dd] = np.mean(omb)
    omb_rmse[dd] = np.sqrt( np.mean( omb**2) )
    omb_rmss[dd] = np.sqrt(np.mean(all_obs_seq_final[dd]['prior ensemble spread']**2 + obs_err_sigma**2) )

    oma =  (all_obs_seq_final[dd]['observations'] - all_obs_seq_final[dd]['posterior ensemble mean'])[flags]
    # print( oma.min(), oma.max()) 
    oma_rmse[dd] = np.sqrt( np.mean( oma**2 ) )
    oma_bias[dd] = np.mean( oma )

    amb = (all_obs_seq_final[dd]['posterior ensemble mean'] - all_obs_seq_final[dd]['prior ensemble mean'])[flags]
    # print( np.mean(amb) )

# --- End of loop over dates


# Plot out metrics
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,2))

ax.axhline(0, color='lightgray')
ax.plot( date_list, omb_bias,'-r')
# ax.plot( date_list, oma_bias,'-k', label='OmA bias')
ax.set_ylabel('6.25um BT Diff (K)')

plt.legend()

ax.xaxis.set_major_formatter(myFmt)
ax.set_title('d) Domain-averaged Innovation of '+prettier_exptname(expt_name), loc='left')



fig.subplots_adjust( left=0.15, bottom=0.2, right=0.97)
plt.savefig('figs/seviri_ch08_%s_metrics.svg' % expt_name)

