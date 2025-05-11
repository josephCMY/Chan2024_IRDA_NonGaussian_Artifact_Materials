'''
Script to plot out how DA influences clouds
-------------------------------------------

1) Does DA cause clouds to have lower RH?
2) Does DA-created clouds have low RH?
3) Does DA-surviving clouds have low RH?
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

from pylib_expt_specific_settings import convert_old_exptname_to_manuscript_name
from copy import deepcopy



''' User control parameters '''
date_st=datetime.datetime.strptime( sys.argv[1], "%Y%m%d%H%M")
date_ed=datetime.datetime.strptime( sys.argv[2], "%Y%m%d%H%M")
exptlist = [sys.argv[3],'NoDA']

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
nZ = 44     # Number of model levels

# Relative humidity bins
rh_bin_edges = np.linspace(0,200,201)/100.
rh_bin_ctrs = (rh_bin_edges[1:]+rh_bin_edges[:-1])/2.
nBins = len(rh_bin_edges)-1




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
      


''' Set up diagnostics container to hold all diagnostics across time and expts'''
diagnostics = {}
for expt in exptlist:
    diagnostics[expt] = {}
    for dname in ['xf cloud rh', 'xa cloud rh', 'xa_rh of DA-created clouds', 'xa_rh of DA-survived clouds', 'xf_rh of DA-survived clouds']:
        diagnostics[expt][dname] = np.zeros( [nD, nZ, nBins] ) + np.nan
    


''' Import diagnostics '''
for expt in exptlist:
    for dd in np.arange(nD):
        date = datelist[dd]
        # Load pickle file containing diagnostics
        pkl_fname = 'pkl_files/%s/%s_CloudyRHStats/diagnostics.pkl' % (expt, date.strftime('%Y%m%d%H%M'))
        # If file does not exist, skip
        if ( not os.path.exists( pkl_fname ) ):
            print( 'missing', pkl_fname)
            continue
        # End of existence check
        pkl_f = open( pkl_fname, 'rb' )
        tmp = pickle.load( pkl_f )
        # print( tmp.keys())
        pkl_f.close()

        # Load diagnostics
        for dname in ['xf cloud rh', 'xa cloud rh', 'xa_rh of DA-created clouds', 'xa_rh of DA-survived clouds', 'xf_rh of DA-survived clouds']:
            diagnostics[expt][dname][dd,:,:] = tmp[dname][:,:]
          


print('Imported RH stats')














'''
    Set up figure
'''
# Init figure
# fig, axes = plt.subplots( nrows = 5, ncols = 1, figsize= (6, 8 ) )
fig = plt.figure( figsize = (6,8) )
axes = [ None, None, None, None, None]
axes[0] = fig.add_axes([0.12, 0.83, 0.93,0.13])
axes[1] = fig.add_axes([0.12, 0.65,0.93,0.13])
axes[2] = fig.add_axes([0.12, 0.42, 0.93,0.13])
axes[3] = fig.add_axes([0.12, 0.24, 0.93,0.13])
axes[4] = fig.add_axes([0.12, 0.06, 0.93,0.13])



# , axes = plt.subplots( nrows = 4, ncols = 1, figsize= (6, 6 ) )
datemesh, pmesh = np.meshgrid( datelist, plvls)




# # Replace top row subplots with one big subplot
# fig.delaxes( axes[0,0])
# fig.delaxes( axes[0,1])
# ax_long = fig.add_axes([0.15,0.7,0.9,0.2])





''' Number of clouds before and after DA as a function of height (left column)'''
expt = exptlist[0]



# DA-created clouds VS DA-survived clouds
ax=axes[0]
prior_clouds = np.sum( diagnostics[expt]['xf cloud rh'], axis = -1 )
poste_clouds = np.sum( diagnostics[expt]['xa cloud rh'], axis = -1 )
diff = (poste_clouds ) / ( prior_clouds+1e-6)  
crange = np.linspace( 0.8,1.6,5)
cnf = ax.contourf( datemesh, pmesh, diff.T, crange, cmap='Greys', extend='max' )
# cnf = ax.contourf( datemesh, pmesh, diff.T, locator=ticker.LogLocator(), cmap='inferno_r', extend='max' )
# ax.contour( datemesh, pmesh, diff.T, [1., 1.1,1.2], colors='k', linestyles=[":","--","-"])
ax.set_title('a) (No. of Analysis Clouds) / (No. of Forecast Clouds)', loc='left')
ax.set_ylim([ 950, 150] )
ax.set_ylabel('Pres. (hPa)')
ax.xaxis.set_major_formatter(myFmt)
cbar = fig.colorbar(cnf, ax=ax)


# DA-created clouds VS DA-survived clouds
ax=axes[1]
da_created_clouds = np.sum( diagnostics[expt]['xa_rh of DA-created clouds'], axis=-1)
da_survived_clouds = np.sum( diagnostics[expt]['xa_rh of DA-survived clouds'], axis=-1)+1e-6
created_frac = (da_created_clouds ) / ( da_survived_clouds + da_created_clouds)  
crange = np.linspace( 0,0.6,7)
cnf = ax.contourf( datemesh, pmesh, created_frac.T, crange, cmap='inferno_r', extend='max' )
# cnf = ax.contourf( datemesh, pmesh, diff.T, locator=ticker.LogLocator(), cmap='inferno_r', extend='max' )
# ax.contour( datemesh, pmesh, diff.T, [1., 1.1,1.2], colors='k', linestyles=[":","--","-"])
ax.set_title('b) (No. of DA-created Clouds) / (No. of Analysis Clouds)', loc='left')
ax.set_ylim([ 950, 150] )
ax.set_ylabel('Pres. (hPa)')
ax.xaxis.set_major_formatter(myFmt)
cbar = fig.colorbar(cnf, ax=ax)



# What is the RH "PDF" like for the prior clouds? 
ax = axes[2]
prior_rh_pdf = np.zeros( [nD, nZ, nBins] )
rh_interval = (rh_bin_edges[1]-rh_bin_edges[0])
for dd in range(nD):
    for kk in range(nZ):
        # print(diagnostics['NoDA']['xf cloud rh'][dd,kk])
        sum_of_clouds = np.sum( diagnostics['NoDA']['xf cloud rh'][dd,kk], axis=-1) +1e-6
        prior_rh_pdf[dd, kk,:] = diagnostics['NoDA']['xf cloud rh'][dd,kk,:] / (sum_of_clouds * rh_interval)
prior_rh_pdf = np.mean(prior_rh_pdf, axis=0)

prior_rh_cdf = np.cumsum( prior_rh_pdf*rh_interval, axis = -1)
# print( prior_rh_cdf.min(), prior_rh_cdf.max() )



prior_cloud_rh_cnf = ax.contourf( rh_bin_ctrs, plvls, prior_rh_cdf, 
                                  np.linspace(0.05,0.95,4),
                                  cmap='viridis_r', extend='max' )
# prior_cloud_rh_cnf = ax.contourf( rh_bin_ctrs, plvls, prior_rh_pdf+1e-6, 
#                                   np.linspace(1,7,7),
#                                   cmap='viridis_r', extend='max')
ax.contour( rh_bin_ctrs, plvls, prior_rh_cdf, np.linspace(0.35,0.95,3), colors='r', linestyles=[':','--','-'], linewidths=2.0)
# ax.contour( rh_bin_ctrs, plvls, prior_rh_cdf, [0.3], colors='r', linestyles='--', linewidths=2)
ax.set_title( 'c) CDF of S in NoDA Clouds', loc='left')
ax.set_ylim([ 950, 150] )
ax.set_xlim([ 0.7, 1.3 ] )
# ax.set_xlabel('Relative Humidity (Pa / Pa)')
ax.set_ylabel('Pres. (hPa)')
cbar = fig.colorbar(prior_cloud_rh_cnf, ax=ax)


# # What is the impact of DA on RH of clouds?
# ax = axes[1]
# poste_rh_pdf = np.zeros( [nD, nZ, nBins] )
# for dd in range(nD):
#     for kk in range(nZ):
#         sum_of_clouds = np.sum( diagnostics[expt]['xa cloud rh'][dd,kk], axis=-1)+1e-6
#         poste_rh_pdf[dd, kk,:] = diagnostics[expt]['xa cloud rh'][dd,kk,:] / (sum_of_clouds * rh_interval)
# poste_rh_pdf = np.mean( poste_rh_pdf, axis=0)
# poste_rh_cdf = np.cumsum( poste_rh_pdf * rh_interval, axis=-1)

# incre_rh_pdf = poste_rh_pdf - prior_rh_pdf
# diff_cloud_cnf = ax.contourf( rh_bin_ctrs, plvls, poste_rh_cdf, np.linspace(0.35,0.95,3), cmap='viridis_r', extend='max')
# ax.contour( rh_bin_ctrs, plvls, prior_rh_cdf, np.linspace(0.35,0.95,3), colors='r', linestyles=[':','--','-'], linewidths=2.0)
# # ax.contour( rh_bin_ctrs, plvls, prior_rh_cdf, [0.3], colors='white', linestyles='-', linewidths=3)
# # ax.contour( rh_bin_ctrs, plvls, prior_rh_cdf, [0.3], colors='r', linestyles='--', linewidths=2)
# ax.set_title( 'b) CDF of RH of Analysis Clouds in %s' % convert_old_exptname_to_manuscript_name(expt), loc='left')
# ax.set_ylim([ 950, 150] )
# ax.set_xlim([ 0.4, 1.05] )
# ax.set_ylabel('Pres. (hPa)')
# # ax.set_xlabel('Relative Humidity (Pa / Pa)')
# fig.colorbar(diff_cloud_cnf, ax=ax)







# What is the RH of DA-created clouds?
ax = axes[3]
fake_cloud_rh_pdf = np.zeros( [nD, nZ, nBins] )
for dd in range(nD):
    for kk in range(nZ):
        sum_of_clouds = np.sum( diagnostics[expt]['xa_rh of DA-created clouds'][dd,kk], axis=-1)+1e-6
        fake_cloud_rh_pdf[dd, kk,:] = diagnostics[expt]['xa_rh of DA-created clouds'][dd,kk,:] / (sum_of_clouds * rh_interval)
fake_cloud_rh_pdf = np.mean( fake_cloud_rh_pdf, axis=0)
fake_cloud_rh_cdf = np.cumsum( fake_cloud_rh_pdf * rh_interval, axis=-1)

fake_cloud_cnf = ax.contourf( rh_bin_ctrs, plvls, fake_cloud_rh_cdf, 
                              np.linspace(0.35,0.95,3), cmap = 'viridis_r', extend='max')
ax.contour( rh_bin_ctrs, plvls, prior_rh_cdf, np.linspace(0.35,0.95,3), colors='r', linestyles=[':','--','-'], linewidths=2.0)
# ax.contour( rh_bin_ctrs, plvls, prior_rh_cdf, [0.3], colors='white', linestyles='-', linewidths=3)
# ax.contour( rh_bin_ctrs, plvls, prior_rh_cdf, [0.3], colors='r', linestyles='--', linewidths=2)
ax.set_title( 'd) CDF of S of DA-Created Analysis Clouds in %s' % convert_old_exptname_to_manuscript_name(expt), loc='left')
ax.set_ylim([ 950, 150] )
ax.set_xlim([ 0.7, 1.3 ] )
ax.set_ylabel('Pres. (hPa)')
# ax.set_xlabel('Relative Humidity (Pa / Pa)')
fig.colorbar(fake_cloud_cnf, ax=ax)




# For DA-surviving clouds, what is the change in RH?
ax = axes[4]
# phys_cloud_prior_rh_pdf = np.zeros( [nD, nZ, nBins] )
# for dd in range(nD):
#     for kk in range(nZ):
#         sum_of_clouds = np.sum( diagnostics[expt]['xf_rh of DA-survived clouds'][dd,kk], axis=-1)+1e-6
#         phys_cloud_prior_rh_pdf[dd, kk,:] = diagnostics[expt]['xf_rh of DA-survived clouds'][dd,kk,:] / (sum_of_clouds * rh_interval)
# phys_cloud_prior_rh_pdf = np.mean( phys_cloud_prior_rh_pdf, axis=0 )

phys_cloud_poste_rh_pdf = np.zeros( [nD, nZ, nBins] )
for dd in range(nD):
    for kk in range(nZ):
        sum_of_clouds = np.sum( diagnostics[expt]['xa_rh of DA-survived clouds'][dd,kk], axis=-1)+1e-6
        phys_cloud_poste_rh_pdf[dd, kk,:] = diagnostics[expt]['xa_rh of DA-survived clouds'][dd,kk,:] / (sum_of_clouds * rh_interval)
phys_cloud_poste_rh_pdf = np.mean( phys_cloud_poste_rh_pdf, axis=0 )
phys_cloud_poste_rh_cdf = np.cumsum( phys_cloud_poste_rh_pdf*rh_interval, axis=-1 )

# change_in_phys_cloud_rh_pdf = phys_cloud_poste_rh_pdf - phys_cloud_prior_rh_pdf
# change_in_phys_cloud_rh_cdf = np.cumsum( change_in_phys_cloud_rh_pdf * rh_interval, axis=-1)
change_in_phys_cloud_rh_pdf_cnf = ax.contourf( rh_bin_ctrs, plvls, phys_cloud_poste_rh_cdf, 
                              np.linspace(0.35,0.95,3), cmap = 'viridis_r', extend='max')
ax.contour( rh_bin_ctrs, plvls, prior_rh_cdf, np.linspace(0.35,0.95,3), colors='r', linestyles=[':','--','-'], linewidths=2.0)
# ax.contour( rh_bin_ctrs, plvls, prior_rh_cdf, [0.3], colors='white', linestyles='-', linewidths=3)
# ax.contour( rh_bin_ctrs, plvls, prior_rh_cdf, [0.3], colors='r', linestyles='--', linewidths=2)
ax.set_title( 'e) CDF of S of DA-Survived Analysis Clouds in %s' % convert_old_exptname_to_manuscript_name(expt), loc='left')
ax.set_ylim([ 950, 150] )
ax.set_xlim([ 0.7, 1.3 ] )
ax.set_xlabel('Saturation Ratio S (Pa / Pa)')
ax.set_ylabel('Pres. (hPa)')
fig.colorbar(change_in_phys_cloud_rh_pdf_cnf, ax=ax)


axes[0].xaxis.set_major_formatter(myBlankFmt)
for ax in axes[2:-1]:
    ax.xaxis.set_major_formatter(myBlankFmt)

# fig.subplots_adjust(hspace=0.5, wspace=0.6, left=0.12, right=1.05, bottom=0.075, top=0.95 )
print( expt, convert_old_exptname_to_manuscript_name(expt))
plt.savefig('figs/CloudInfo_%s.pdf' % convert_old_exptname_to_manuscript_name(expt))
plt.close()


# Extra
crange = np.linspace( 0.5,1.0,11)
print( prior_clouds.shape)
cnf = plt.contourf( datemesh, pmesh, (da_survived_clouds/prior_clouds).T, crange, cmap='inferno_r' )
plt.colorbar(cnf)
plt.ylim([950,150])
plt.title('(No. of DA-survived Clouds) / (No. of Forecast Clouds)')
plt.savefig('figs/da_survived_fraction_%s.pdf' % convert_old_exptname_to_manuscript_name(expt))