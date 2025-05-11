'''
Script to plot out validation statistics during cycling.
--------------------------------------------------------
Plots to make:
    1) NoDA-relative RMSEs as a function of level and time.
    2) Normalized biases as a function of level and time
'''

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import math
import datetime
import sys
import pickle
from netCDF4 import Dataset as ncopen
import matplotlib.dates as mdates
from matplotlib import ticker
myFmt = mdates.DateFormatter('%b-%d\n%HUTC')
myBlankFmt = mdates.DateFormatter(' ')
import os
from scipy.stats import spearmanr
from pylib_expt_specific_settings import *
from copy import deepcopy
from math import pi as PI


''' User control parameters '''
date_st=datetime.datetime.strptime( '202301130000', "%Y%m%d%H%M")
date_ed=datetime.datetime.strptime( '202301160000', "%Y%m%d%H%M")
exptlist = ['FCSTCLOUD']
#['BASIC', 'HIPASS','NoCloudUpdate','HIPASS_NoCloudUpdate'] 

vnamelist = ['T', 'QVAPOR' ] #,'QVAPOR','P','PH'] #,'QRAIN','QSNOW','QICE','QGRAUP']
name_list = {}
for expt in exptlist:
  name_list[expt] = expt


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


stagelist = ['xf','xa']
t_int = datetime.timedelta( minutes=60 )
nK = 201    # Number of wavenumbers
nZ = 44 # Number of model levels

## Special handling for PH
#if vname == 'PH':
#  nZ = 45

''' Add on NoDA as an expt '''
exptlist.insert(0,'NoDA')


''' Set up list of dates '''
datelist = []
date = date_st + datetime.timedelta(minutes=0)
while date <= date_ed:
    datelist.append(date)
    date += t_int
datelist = np.array(datelist)
nD = len(datelist)


''' Figuring out dates to label '''
dateticks = [ datelist[0], datelist[0]+ (datelist[-1]-datelist[0])/3, 
              datelist[0]+ (datelist[-1]-datelist[0])*2/3, \
              datelist[-1] ]
      












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
    LOAD OSSE DART OBS DATA
'''
# Init array to hold all data entries
all_obs_seq_final = [None for i in range(nD)]

expt = exptlist[1]

# Loop over all dates
for dd in range(nD):

    date_nw = datelist[dd]

    # Load real data 
    fname = date_nw.strftime('DART_files/'+expt+'/obs_seq.final.%Y%m%d%H%M')
    f = open( fname, 'r' )
    all_obs_seq_final[dd] = read_obs_seq_final_obs_and_prior( f, obs_kind=260 )
    f.close()

# --- End of loop over all dates








'''
    COMPUTE INNOVATION BIASES
'''
omb_bias = np.zeros(nD)
# Compute metrics for all dates
for dd in range(nD):

    # Removing strange values...
    flags = (all_obs_seq_final[dd]['posterior ensemble mean'] < 100) + (all_obs_seq_final[dd]['posterior ensemble mean'] > 350 )
    flags += (all_obs_seq_final[dd]['observations'] < 100) + (all_obs_seq_final[dd]['observations'] > 350 )
    flags = (flags == 0)

    omb = (all_obs_seq_final[dd]['observations'] - all_obs_seq_final[dd]['prior ensemble mean'])[flags]
    omb_bias[dd] = np.mean(omb)

# -- End of loop




















''' Set up diagnostics container to hold all diagnostics across time, channel and expts'''
diagnostics = {}
for vn in vnamelist:
  diagnostics[vn] = {}
  for expt in exptlist:
    diagnostics[vn][expt] = {}
    for ss in stagelist:
      diagnostics[vn][expt][ss] = {}
      for dname in ['EmT_bias','MSS','MSE','Dom Avg']:
        diagnostics[vn][expt][ss][dname] = np.zeros( [nD, nZ], dtype = float ) + np.nan
    diagnostics[vn][expt]['xt Dom Avg'] = np.zeros( [nD, nZ], dtype = float ) + np.nan


''' Import diagnostics '''
for vn in vnamelist:
  for expt in exptlist:
    for dd in np.arange(nD):
      date = datelist[dd]
      # Load pickle file containing diagnostics
      pkl_fname = 'pkl_files/%s/%s_%s/validation.pkl' % (expt, date.strftime('%Y%m%d%H%M'), vn )
      # If file does not exist, skip
      if ( not os.path.exists( pkl_fname ) ):
        print( 'missing', pkl_fname)
        continue
      # End of existence check
      pkl_f = open( pkl_fname, 'rb' )
      tmp = pickle.load( pkl_f )
      pkl_f.close()
      # Load mean field
      diagnostics[vn][expt]['xt Dom Avg'][dd] = tmp['xt Dom Avg']
      # Load stage-specific information
      for dname in ['EmT_bias','MSS','MSE','Dom Avg']:
        for ss in stagelist:
          if vn == 'PH' or vn == 'W':
            diagnostics[vn][expt][ss][dname][dd] = tmp[ss][dname][1:]
          else:
            diagnostics[vn][expt][ss][dname][dd] = tmp[ss][dname]


# Special handling for mixing ratios to convert to g/kg
for vn in vnamelist:
  if vn[0] == 'Q':
    for expt in exptlist:
      diagnostics[vn][expt]['xt Dom Avg'] *= 1e3
      for ss in stagelist:
        diagnostics[vn][expt][ss]['EmT_bias'] *= 1e3
        diagnostics[vn][expt][ss]['Dom Avg'] *= 1e3
        for dname in ['MSS','MSE']:
          diagnostics[vn][expt][ss][dname] *= 1e6

print('Imported RMS stats')














































letters = ['a','b','c','d']


''' Evolution of biases '''

# Plotting out biases
datemesh, pmesh = np.meshgrid( datelist[1:], plvls)

# Init figure
fig, axes = plt.subplots( nrows = 3, ncols = 2, figsize= (6, 6 )  )

# Adjust subplots
fig.subplots_adjust( top = 0.9, bottom = 0.35, left=0.125, right=0.925,
                      hspace = 0.4 , wspace=0.2 )


# Add on row for innovation statistics
ax_inno = fig.add_subplot([0.125,0.08,0.8,0.15])




# Color ranges
crange = {'T': np.linspace(-0.05,0.05,6), "QVAPOR": np.linspace(-0.02,0.02,6) }


# Loop over various variables
for vv in range(len(vnamelist)):

  # Variable name
  vn = vnamelist[vv]

  # Compute impact of dom-avg analysis incre on bias at time t
  incre_impact_on_bias = diagnostics[vn][expt]['xa']['Dom Avg'][:-1] - diagnostics[vn][expt]['xf']['Dom Avg'][:-1]

  # Compute domain-averaged model response impact on bias
  # 1) Compute change in bias from analysis at time t to forecast at time t+tau
  model_impact_on_bias = diagnostics[vn][expt]['xf']['Dom Avg'][1:] - diagnostics[vn][expt]['xa']['Dom Avg'][:-1] 
  # 2) Account for evolution of truth state in time
  model_impact_on_bias -= ( diagnostics[vn][expt]['xt Dom Avg'][1:] -  diagnostics[vn][expt]['xt Dom Avg'][:-1] )

  # Normalizer to enable multi-layer plot
  normalizer = np.sqrt( np.mean(diagnostics[vn]['NoDA']['xf']['MSE'][:], axis=0) )


  # Domain-averaged bias change from forecast at t to forecast at t+tau
  overall_bias_change = diagnostics[vn][expt]['xf']['Dom Avg'][1:] - diagnostics[vn][expt]['xf']['Dom Avg'][:-1]
  overall_bias_change -= ( diagnostics[vn][expt]['xt Dom Avg'][1:] - diagnostics[vn][expt]['xt Dom Avg'][:-1] )

  check = np.abs(overall_bias_change - (incre_impact_on_bias + model_impact_on_bias) )/normalizer
  print( check.min(),check.max() )
  # print( (diagnostics[vn][expt]['xt Dom Avg'][1:] - diagnostics[vn][expt]['xt Dom Avg'][:-1]) /)
  



  # Plot domain-averaged analysis increment impact on bias
  bias_impact = incre_impact_on_bias / normalizer
  cnf1 = axes[0,vv].contourf( datemesh, pmesh, bias_impact.T, crange[vn], cmap = 'RdBu_r', extend='both' )

  # Formatting for analysis increment                
  axes[0,vv].set_title('a' + str(vv+1) +') ' +vn+' Incre. $\Delta$sBias',  loc='left')
  plt.colorbar(cnf1, ax = axes[0,vv])



  # Plot domain-averaged model response impact on bias
  bias_impact = model_impact_on_bias / normalizer
  cnf1 = axes[1,vv].contourf( datemesh, pmesh, bias_impact.T, crange[vn], cmap = 'RdBu_r', extend='both' )

  # Formatting for analysis increment                      
  axes[1,vv].set_title('b' + str(vv+1) +') ' +vn+' Resp. $\Delta$sBias',  loc='left')
  plt.colorbar(cnf1, ax = axes[1,vv])



  # Plot overall evolution of bias
  overall_bias_change = overall_bias_change / normalizer
  cnf1 = axes[2,vv].contourf( datemesh, pmesh, overall_bias_change.T, crange[vn], cmap = 'RdBu_r', extend='both' )

  # Formatting for analysis increment                      
  axes[2,vv].set_title('c' + str(vv+1) +') ' +vn+' $\Delta$sBias',  loc='left')
  plt.colorbar(cnf1, ax = axes[2,vv])


  # Cosmetic stuff
  for i in range(3):
    axes[i,vv].set_xticks(dateticks)
    axes[i,vv].set_ylim([ 1000, 20.] )
    axes[i,vv].set_yticks([1000,750,500,250])
    axes[i,0].set_ylabel('Pres. (hPa)')
    axes[i,1].set_ylabel(' ')
    axes[i,1].yaxis.set_major_formatter(myBlankFmt)

    # Handling axis ticks
    if i < axes.shape[0]-1:
      axes[i,vv].xaxis.set_major_formatter(myBlankFmt)
    else:
      axes[i,vv].xaxis.set_major_formatter(myFmt)




# Throw on subplot about innovation biases
ax_inno.axhline(0, color='lightgray')
ax_inno.plot( datelist, omb_bias,'-r')
ax_inno.set_ylabel('Inno. (K)')
ax_inno.xaxis.set_major_formatter(myFmt)
ax_inno.set_title('d) Domain-averaged Innovation of '+convert_old_exptname_to_manuscript_name(expt), loc='left')




# Label figuew, v
plt.suptitle( 'Evolution of Biases for %s' % ( convert_old_exptname_to_manuscript_name(expt)) ) 
  

# Save figure
plt.savefig('figs/bias_evolution_%s.pdf' % (convert_old_exptname_to_manuscript_name(expt)))
plt.close()