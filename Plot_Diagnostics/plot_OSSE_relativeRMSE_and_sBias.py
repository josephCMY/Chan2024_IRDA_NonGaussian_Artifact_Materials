'''
Script to plot out validation statistics during cycling.
--------------------------------------------------------
Plots to make:
    1) Relative RMSEs as a function of level and time.
    2) Scaled biases as a function of level and time
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




''' User control parameters '''
# Start date
date_st=datetime.datetime.strptime( sys.argv[1], "%Y%m%d%H%M")

# End date
date_ed=datetime.datetime.strptime( sys.argv[2], "%Y%m%d%H%M")

# Experiments to compare. 
exptlist = [sys.argv[3], sys.argv[4]]



vnamelist = ['U','V','T', 'QVAPOR' ]
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



# # Normalize MSE
# for vn in vnamelist:
#   for expt in exptlist:

#     # Normalizing!
#     for ss in stagelist:
#       diagnostics[vn][expt][ss]['norm RMSE']  \
#       = ( np.sqrt( diagnostics[vn][expt][ss]['MSE'] )
#           / np.sqrt( diagnostics[vn][noda_name][ss]['MSE'] )
#         )
#     # End of loop over stages
#   # End of loop over expt names
# # End of loop over variable names





''' Normalize biases using NoDA data '''
# Normalize biases 
for vn in vnamelist:
  for expt in exptlist:

    # Normalizing!
    for ss in stagelist:
      diagnostics[vn][expt][ss]['norm bias']  \
      = ( diagnostics[vn][expt][ss]['EmT_bias'] 
           / np.mean( np.sqrt( diagnostics[vn][noda_name][ss]['MSE'] ), axis=0 )
        )
    # End of loop over stages
  # End of loop over expt names
# End of loop over variable names





letters = ['a','b','c','d']


''' Plot RMSEs as a function of cycling time and model level '''
# Plotting out RMSD and all for the chosen wavebands
datemesh, pmesh = np.meshgrid( datelist, plvls)
for ss in stagelist:
  # Init figure
  fig, axes_all = plt.subplots( nrows = len(vnamelist), ncols = 2, figsize= (5, 6 )  )

  # Plot out RMSEs
  axes = axes_all[:,0]
  for vv in range(len(vnamelist)):
    vn =vnamelist[vv]
    expt1 = exptlist[-2]
    expt2 = exptlist[-1]
    # NoDA-relative RMSes
    cnf1 = axes[vv].contourf( datemesh, pmesh,
                              np.sqrt(diagnostics[vn][expt2][ss]['MSE'].T / diagnostics[vn][expt1][ss]['MSE'].T), 
                              # ( diagnostics[vn][expt][ss]['norm RMSE'] - diagnostics[vn]['NoDA'][ss]['norm RMSE'] ).T,
                              np.linspace( 0.7, 1.3, 6),
                              cmap = 'PRGn_r', extend='both',
                              zorder=0
                            )


    # Consistency ratios
    cr = np.sqrt( diagnostics[vn][expt2][ss]['MSS']
                  / diagnostics[vn][expt2][ss]['MSE']
                  ).T
    # cn = axes[ee].contour( datemesh, pmesh, cr,
    #                         np.arange(6)*0.2, 
    #                         colors=    ['r','r','r','k','k','k'],
    #                         linestyles=[':','--','-',':','--','-']
    # )
    cn1 = axes[vv].contourf( datemesh, pmesh, cr,
                            [0.5,0.75], colors = 'none',
                            hatches = ['x'] )
    cn2 = axes[vv].contourf( datemesh, pmesh, cr,
                            [0.,0.5], colors = 'none',
                            hatches = ['xxx'] )
    for collection in cn1.collections:
        collection.set_edgecolor( 'r' )
    for collection in cn2.collections:
        collection.set_edgecolor( 'r' )

    # axes[vv].contour( datemesh, pmesh, cr,
    #                     [0.5,0.75], colors='r')


                      
    axes[vv].set_title(letters[vv]+'1) '+vn+' rel. RMSE', loc='left')
    axes[vv].set_xticks(dateticks)
    axes[vv].set_ylim([ 1000, 20.] )
    axes[vv].set_yticks([1000,750,500,250])
    axes[vv].set_ylabel('Pres. (hPa)')
    plt.colorbar(cnf1,ax=axes[vv])

    
    # Handling axis ticks
    if vv < len(vnamelist)-1:
      axes[vv].xaxis.set_major_formatter(myBlankFmt)
    else:
      axes[vv].xaxis.set_major_formatter(myFmt)




  # Plot out biases
  axes = axes_all[:,1]
  for vv in range(len(vnamelist)):
    vn =vnamelist[vv]
    expt = exptlist[-1]
    # NoDA-relative biases
    cnf1 = axes[vv].contourf( datemesh, pmesh,
                              ( diagnostics[vn][expt][ss]['norm bias'] ).T,
                              np.linspace( -0.5,0.5, 6),
                              cmap = 'RdBu_r', extend='both',
                              zorder=0
                            ) 

                      
    axes[vv].set_title(letters[vv]+'2) '+vn+' sBias',  loc='left')
    axes[vv].set_xticks(dateticks)
    axes[vv].set_ylim([ 1000, 20.] )
    axes[vv].set_yticks([1000,750,500,250])
    plt.colorbar(cnf1,ax=axes[vv])
    axes[vv].set_ylabel('')
    axes[vv].yaxis.set_major_formatter(myBlankFmt)

    # Handling axis ticks
    if vv < len(vnamelist)-1:
      axes[vv].xaxis.set_major_formatter(myBlankFmt)
    else:
      axes[vv].xaxis.set_major_formatter(myFmt)




  # Adjust subplots
  fig.subplots_adjust( top = 0.85, bottom = 0.07, left=0.13, right=0.95,
                        hspace = 0.4 , wspace=0.2)
  # Label figuew, v
  plt.suptitle( "%s's %s-relative RMSEs (left)\n and %s's sBiases (right)"  
        % ( 
            exptlist[-1],
            exptlist[-2],
            exptlist[-1]
          ),
          linespacing=1.5,
          fontsize=11.5
      )
  

    

  # Save figure
  #plt.savefig('figs/rmse_%s_%s.png' % (ss, vn), dpi=500, bbox_to_inches='tight')
  plt.savefig('figs/norm_metrics_%s_vs_%s.pdf' % (exptlist[-2], exptlist[-1]))
  plt.close()




