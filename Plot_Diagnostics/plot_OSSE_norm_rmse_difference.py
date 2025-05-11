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
myFmt = mdates.DateFormatter('%b-%d\n%HUTC')
myBlankFmt = mdates.DateFormatter(' ')
import os
from scipy.stats import spearmanr
from pylib_expt_specific_settings import *



''' User control parameters '''
date_st=datetime.datetime.strptime( sys.argv[1], "%Y%m%d%H%M")
date_ed=datetime.datetime.strptime( sys.argv[2], "%Y%m%d%H%M")
exptlist = [sys.argv[3], sys.argv[4]]
#['BASIC', 'HIPASS','NoCloudUpdate','HIPASS_NoCloudUpdate'] 
expt0 = exptlist[0]
expt1 = exptlist[1]

vnamelist = ['U','V','T', 'QVAPOR' ] #,'QVAPOR','P','PH'] #,'QRAIN','QSNOW','QICE','QGRAUP']
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


stagelist = ['xf']
t_int = datetime.timedelta( minutes=60 )
nK = 201    # Number of wavenumbers
nZ = 44 # Number of model levels

## Special handling for PH
#if vname == 'PH':
#  nZ = 45

''' Add on NoDA as an expt '''
if 'NoDA' not in exptlist:
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
      


''' Set up diagnostics container to hold all diagnostics across time, channel and expts'''
diagnostics = {}
for vn in vnamelist:
  diagnostics[vn] = {}
  for expt in exptlist:
    diagnostics[vn][expt] = {}
    for ss in stagelist:
      diagnostics[vn][expt][ss] = {}
      for dname in ['EmT_bias','MSS','MSE']:
        diagnostics[vn][expt][ss][dname] = np.zeros( [nD, nZ], dtype = float ) + np.nan


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
      for dname in ['EmT_bias','MSS','MSE']:
        for ss in stagelist:
          if vn == 'PH' or vn == 'W':
            diagnostics[vn][expt][ss][dname][dd] = tmp[ss][dname][1:]
          else:
            diagnostics[vn][expt][ss][dname][dd] = tmp[ss][dname]


# Special handling for mixing ratios to convert to g/kg
for vn in vnamelist:
  if vn[0] == 'Q':
    for expt in exptlist:
      for ss in stagelist:
        diagnostics[vn][expt][ss]['EmT_bias'] *= 1e3
        for dname in ['MSS','MSE']:
          diagnostics[vn][expt][ss][dname] *= 1e6

print('Imported RMS stats')




''' Normalize MSEs using NoDA data '''
# Seek NoDA expt name
for expt in exptlist:
  if ( 'NoDA' in expt ):
    noda_name = expt
    break
  # endif
# end seeking loop


# Compute norm factors
norm_factor = {}
for vn in vnamelist:
  norm_factor[vn] = {}
  for ss in stagelist:
    norm_factor[vn][ss] = np.mean( 
        np.sqrt( diagnostics[vn][noda_name][ss]['MSE'] ),
        axis=0
    )


# Normalize MSE
for vn in vnamelist:
  for expt in exptlist:
    # Normalizing!
    for ss in stagelist:
      diagnostics[vn][expt][ss]['norm RMSE'] = (
          np.sqrt( diagnostics[vn][expt][ss]['MSE'] ) 
          / np.sqrt( diagnostics[vn][noda_name][ss]['MSE'] )
      )
    # End of loop over stages
  # End of loop over expt names
# End of loop over variable names




''' Normalize biases using NoDA data '''
# Normalize biases 
for vn in vnamelist:
  for expt in exptlist:

    # Normalizing!
    for ss in stagelist:
      diagnostics[vn][expt][ss]['norm bias']  \
      = ( diagnostics[vn][expt][ss]['EmT_bias'] 
           / norm_factor[vn][ss]
          )
    # End of loop over stages
  # End of loop over expt names
# End of loop over variable names







letters = ['a','b','c','d']


''' Plot RMSEs and differences in RMSEs as a function of cycling time and model level '''
# Plotting out RMSD and all for the chosen wavebands
datemesh, pmesh = np.meshgrid( datelist, plvls)
for ss in stagelist:
  # Init figure
  fig, axes_all = plt.subplots( nrows = 2, ncols = 2, figsize= (5, 3)  )
  axes = [val for sub in axes_all for val in sub]


  # Plot out NoDA-relative RMSE 
  for vv in range(len(vnamelist)):
    vn =vnamelist[vv]
    expt = exptlist[1]

    # Relative diff in RMSEs
    diff = ( diagnostics[vn][expt1][ss]['norm RMSE'] ).T / ( diagnostics[vn][expt0][ss]['norm RMSE'] ).T  
    # diff[ np.abs(diff)<0.1 ] = np.nan
    cnf3 = axes[vv].contourf( datemesh, pmesh,
                              diff,
                              np.linspace( 0.5,1.5, 6),
                              cmap = 'PRGn_r', extend='both',
                              zorder=0
                            )
    plt.colorbar( cnf3, ax = axes[vv])


                      
    axes[vv].set_title('%s) %s rel. RMSE ' % 
                       (  letters[vv], vn), 
                    loc='left')
    axes[vv].set_xticks(dateticks)
    axes[vv].set_ylim([ 1000, 20.] )
    axes[vv].set_yticks([1000,750,500,250])
    axes[vv].set_ylabel('Pres. (hPa)')
  
    # # # Blanking out pressure ticks
    # # axes[vv].yaxis.set_major_formatter(myBlankFmt)

  # Only apply date ticks to bottom row
  for i in range(2):
    axes_all[0,i].xaxis.set_major_formatter(myBlankFmt)
    axes_all[1,i].xaxis.set_major_formatter(myFmt)



  # Only applying pressure coordinates ticks to the left column
  for i in range(2):
    axes_all[i,1].yaxis.set_major_formatter(myBlankFmt)
    axes_all[i,1].set_ylabel('')


  # Adjust subplots
  fig.subplots_adjust( top = 0.8, bottom = 0.13, left=0.14, right=0.95,
                        hspace = 0.4 , wspace=0.2)








  plt.suptitle( '%s RMSE divided by %s RMSE ' % (convert_old_exptname_to_manuscript_name(expt1), convert_old_exptname_to_manuscript_name(expt0)) )
  # Label figuew, v
  # plt.suptitle( '%s $\Delta$sRMSEs (left)  and  %s minus %s $\Delta$sRMSEs (right)' % ( convert_old_exptname_to_manuscript_name(expt1), convert_old_exptname_to_manuscript_name(expt1) , convert_old_exptname_to_manuscript_name(expt0)) ) 
    

  # Save figure
  #plt.savefig('figs/rmse_%s_%s.png' % (ss, vn), dpi=500, bbox_to_inches='tight')
  plt.savefig('figs/diff_norm_metrics_%s_%s_%s.pdf' % (ss, convert_old_exptname_to_manuscript_name(expt0), convert_old_exptname_to_manuscript_name(expt1)))
  plt.close()




