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
from pylib_read_obs_priors import *

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
    LOAD USER INPUTS
'''
date_st = datetime.strptime('202301130000', '%Y%m%d%H%M')
date_ed = datetime.strptime('202301160000', '%Y%m%d%H%M')
time_interval = float( 60 )












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
    all_real_obs[dd] = np.array( read_obs_seq( f, obs_kind=260, entryname='observation' )['observation'] )
    f.close()

# --- End of loop over all dates












'''
    COMPUTE OBSERVATION CONSISTENCY RATIOS
'''
expt_name='NoDA'
ens_size=50
obs_err_sigma=1.0


# Dictionary to hold onto CR values
cr_dict = {}
for vname in ['NoDA vs Real Obs', 'NoDA vs Synth Obs','Hypothetical Gaussian']:
    cr_dict[vname] = []

# Dictionary to control plotting styles
plot_dict = {
    'NoDA vs Real Obs': {'lw': 2, 'c':'r', 'zo': 10},
    'NoDA vs Synth Obs': {'lw': 2, 'c':'k', 'zo': 8},
    'Hypothetical Gaussian': {'lw': 2, 'c':'dodgerblue', 'zo': 9}#,
    # 'Std Norm': {'lw': 4, s'c':'darkgray', 'zo': 8},
}



# Loop over all dates
for dd in range(nDates):

    date_nw = date_list[dd]

    # Load NoDA ensemble and synthetic obs
    fname = date_nw.strftime('DART_files/'+expt_name+'/obs_prior/obs_seq.final.%Y%m%d%H%M')
    f = open( fname, 'r' )
    obs_dict = read_obs_seq_prior_and_obsval( f, obs_kind=260, ens_size=ens_size )
    f.close()
    all_obs_prior_names = ['prior ensemble member     %2d' % (ee+1) \
                            for ee in range(ens_size)]
    all_obs_priors = np.array(
        [ obs_dict[all_obs_prior_names[ee]] for ee in range(ens_size) ] 
    )
    synth_obs = np.array( obs_dict['observations' ])

    # Compute signed consistency ratio for every real obs & store in NoDA
    inno = all_real_obs[dd] - np.mean( all_obs_priors, axis=0)
    sigma2 = (
        np.var( all_obs_priors, axis=0, ddof=1)
        + obs_err_sigma**2
    )
    cr_dict['NoDA vs Real Obs'].append(
        np.sign(inno) * np.sqrt(
            inno**2 / sigma2
        )
    )

    # Compute signed consistency ratio for every synthetic obs & store in NoDA
    inno = synth_obs - np.mean( all_obs_priors, axis=0)
    cr_dict['NoDA vs Synth Obs'].append(
        np.sign(inno) * np.sqrt(
            inno**2 / sigma2
        )
    )


    # Generate Gaussian situation (uses NoDA's innovation variance)
    # Assumes unbiased observations
    noda_prior_sigma = np.std( all_obs_priors, axis=0, ddof=1)
    noda_prior_mean = np.mean(all_obs_priors, axis=0)
    gauss_prior_samples = np.random.normal( size = all_obs_priors.shape ) * noda_prior_sigma  + noda_prior_mean 
    gauss_obs_samples = np.random.normal( size = noda_prior_sigma.shape ) * noda_prior_sigma  + noda_prior_mean 
    gauss_obs_samples += np.random.normal( size = noda_prior_sigma.shape ) * obs_err_sigma
    gauss_prior_var = np.var( gauss_prior_samples, axis = 0, ddof=1)
    gauss_prior_mean = np.mean( gauss_prior_samples, axis=0)
    inno = gauss_obs_samples - gauss_prior_mean
    sigma2 = gauss_prior_var + obs_err_sigma**2
    cr_dict['Hypothetical Gaussian'].append(
        np.sign(inno) * np.sqrt(
            (inno**2)/sigma2
        )
    )




# Plot histogram of all signed CRs
fig = plt.figure(figsize=(6,3))
# cr_bin_edges = np.linspace(-50,50,101)
# bin_ctrs = (cr_bin_edges[1:]+cr_bin_edges[:-1])/2

for iv, vname in enumerate(cr_dict.keys()):
    # Plot
    signed_cr = ( np.array( cr_dict[vname] ) ).flatten()
    plt.hist( 
        signed_cr, bins = 101, histtype='step', label=vname, 
        color=plot_dict[vname]['c'], linewidth=plot_dict[vname]['lw'],
        zorder = plot_dict[vname]['zo']
    )

plt.axvline( 0, color='k', linestyle='--')
plt.legend()
plt.ylabel('Occurrence Frequency')
plt.yscale('log')
# plt.xlim([-0.02, 0.3])
plt.xlabel('Sign(d) * CR')
plt.title('Consistency Ratios implied by IR-BT Innovation Distributions')

plt.tight_layout()

plt.savefig('figs/signed_crs.pdf')

# --- End of loop over all dates
