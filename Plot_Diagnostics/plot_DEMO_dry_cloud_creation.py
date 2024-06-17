'''
    CONCEPTUAL DEMONSTRATION TO SHOW HOW DRY CLOUDS CAN BE CREATED BY ENKF
'''


import numpy as np
from matplotlib import use as mpl_use
mpl_use('agg')
import matplotlib.pyplot as plt
from scipy.stats import multivariate_normal



# EnSRF function (mathematically identical to EAKF)
# This function has been sanity checked.
def ensrf( xf_ens, yf_ens, obs, obs_err_var ):

    # Compute kalman gain
    xf_pert = xf_ens - np.mean( xf_ens )
    yf_pert = yf_ens - np.mean( yf_ens )
    ens_size = len(xf_ens)
    cov_xy = np.sum( xf_pert * yf_pert ) / (ens_size - 1.)
    cov_yy = np.sum( yf_pert * yf_pert ) / (ens_size - 1.)
    kalman = cov_xy / (cov_yy + obs_err_var)

    # Compute posterior ensemble
    xa_avg = np.mean(xf_ens) + kalman * ( obs - np.mean(yf_ens) )
    phi = 1 + np.sqrt(  obs_err_var / (cov_yy+obs_err_var)  )
    phi = 1./phi
    xa_pert = xf_pert - kalman * phi * yf_pert
    xa_ens = xa_avg + xa_pert

    return xa_ens, cov_xy/cov_yy




















'''
    Prepare and run demonstration
'''


# Manually prepared prior ensemble
irbt_xf_ens = np.array([
    240, 238, 236,
    225, 215, 210
])

qvapor_xf_ens = np.array([
    0.0016, 0.0018, 0.002,
    0.0025, 0.0025,0.0025
])*1e3

qhydro_xf_ens = np.array([
    0.0, 0.0, 0.0,
    0.0001, 0.0002,0.0005
])*1e3

obs = 220
obs_err_var = (2**2) * 1.




# Construct posteriors
irbt_xa_ens, reg = ensrf( irbt_xf_ens*1., irbt_xf_ens*1., obs, obs_err_var )
qvapor_xa_ens, qvapor_ir_reg = ensrf( qvapor_xf_ens*1., irbt_xf_ens*1., obs, obs_err_var )
qhydro_xa_ens, qhydro_ir_reg = ensrf( qhydro_xf_ens*1., irbt_xf_ens*1., obs, obs_err_var )








'''
    Fit all data into one convenient dictionary
'''
all_data_dict = {}

all_data_dict['ir'] = {}
all_data_dict['ir']['xf'] = irbt_xf_ens
all_data_dict['ir']['xa'] = irbt_xa_ens
all_data_dict['ir']['name'] = 'IR BT'
all_data_dict['ir']['unit'] = 'K'
all_data_dict['ir']['range'] = np.linspace(205,245,101)

all_data_dict['qv'] = {}
all_data_dict['qv']['xf'] = qvapor_xf_ens
all_data_dict['qv']['xa'] = qvapor_xa_ens
all_data_dict['qv']['name'] = 'QVAPOR'
all_data_dict['qv']['unit'] = 'g/kg'
all_data_dict['qv']['range'] = np.linspace(1.5,2.7,101)

all_data_dict['qh'] = {}
all_data_dict['qh']['xf'] = qhydro_xf_ens
all_data_dict['qh']['xa'] = qhydro_xa_ens
all_data_dict['qh']['name'] = 'QHYDROS'
all_data_dict['qh']['unit'] = 'g/kg'
all_data_dict['qh']['range'] = np.linspace(-0.05,0.55,101)










'''
    Plot demonstration

    Row 1: bivariate plots of prior ens, prior Gaussian PDFs and obs likelihood function
    Row 2: bivariate plots of posterior ens & posterior gaussian PDFs
'''
ens_size = len(qhydro_xa_ens)

markers = [r'$\text{A}$',r'$\text{B}$',r'$\text{C}$',r'$\text{D}$',r'$\text{E}$',r'$\text{F}$']
alph = ['a','b','c']

plot_combos = ([
    ['qv', 'ir'],
    ['qh', 'ir'],
    ['qv','qh']
])

fig, axs = plt.subplots( nrows=2, ncols=3, figsize=(8,6.5))

# axs = axs.T



# Plot various combo's fcst rows
for icol, combo in enumerate( plot_combos ):

    # extract variable keys
    v1, v2 = combo

    # Plot forecast distribution
    mean = np.array( [ np.mean(all_data_dict[v1]['xf']), np.mean(all_data_dict[v2]['xf']) ] )
    cov = np.cov( all_data_dict[v1]['xf'], all_data_dict[v2]['xf'], ddof=1)
    rv_obj = multivariate_normal(mean = mean, cov = cov)
    v1mesh, v2mesh = np.meshgrid( all_data_dict[v1]['range'], all_data_dict[v2]['range'])
    pos = np.dstack( (v1mesh, v2mesh) )
    pdf = rv_obj.pdf(pos) 
    pdf /= pdf.max()
    CS = axs[0,icol].contour( v1mesh, v2mesh, pdf, [0.5, 0.9], linestyles=[':','-'], linewidths=1, colors='darkgray', zorder=0)
    # axs[0,icol].clabel(CS, CS.levels, inline=True)
    # cnf1 = axs[0,icol].pcolormesh( v1mesh, v2mesh, pdf, vmin=0.5, vmax=1.2, cmap ='Greys', alpha=0.5, zorder=0 )
    # fig.colorbar( cnf, ax=axs[0,icol])

    # Plot individual members for first row
    for imem in range( ens_size ):
        axs[0,icol].plot( 
            all_data_dict[v1]['xf'][imem], all_data_dict[v2]['xf'][imem], markersize=7, mew = 0.5,
            marker=markers[imem], linewidth=0, color='k', label='Forecast Mem %s' % (markers[imem]) 
        )

    # Label plots
    axs[0,icol].set_xlabel( '%s (%s)' % (all_data_dict[v1]['name'], all_data_dict[v1]['unit']) )
    axs[0,icol].set_ylabel( '%s (%s)' % (all_data_dict[v2]['name'], all_data_dict[v2]['unit']) )
    axs[0,icol].set_title( 'a%d) Forecast %s & %s' % (icol+1, v2.upper(), v1.upper() ), loc='left')


    # Throw on saturation lines
    if v1 == 'qv':
        axs[0,icol].axvline( 2.5, color = 'k', zorder=1, linestyle="--", linewidth=0.5, label='RH = 100%')





# Plot various combo's analysis rows
for icol, combo in enumerate( plot_combos ):

    # extract variable keys
    v1, v2 = combo

    # Plot analysis distribution
    mean = np.array( [ np.mean(all_data_dict[v1]['xa']), np.mean(all_data_dict[v2]['xa']) ] )
    cov = np.cov( all_data_dict[v1]['xa'], all_data_dict[v2]['xa'], ddof=1)
    rv_obj = multivariate_normal(mean = mean, cov = cov)
    v1mesh, v2mesh = np.meshgrid( all_data_dict[v1]['range'], all_data_dict[v2]['range'])
    pos = np.dstack( (v1mesh, v2mesh) )
    pdf = rv_obj.pdf(pos) 
    pdf /= pdf.max()
    CS = axs[1,icol].contour( v1mesh, v2mesh, pdf, [0.5, 0.9], linestyles=[':','-'], linewidths=1, colors='lightcoral', zorder=0)
    # axs[1,icol].clabel(CS, CS.levels, inline=True)
    # cnf2 = axs[1,icol].pcolormesh( v1mesh, v2mesh, pdf, vmin=0.5, vmax =1.2, cmap = 'Reds', alpha=0.5, zorder=0 )
    # fig.colorbar( cnf, ax=axs[1,icol])


    # Plot individual members for first row
    for imem in range( ens_size ):
        axs[1,icol].plot( 
            all_data_dict[v1]['xf'][imem], all_data_dict[v2]['xf'][imem], markersize=7, mew = 0.5,
            marker=markers[imem], linewidth=0, color='k', label='Forecast Mem %s' % (markers[imem]) 
        )
    for imem in range( ens_size ):
        axs[1,icol].plot( 
            all_data_dict[v1]['xa'][imem], all_data_dict[v2]['xa'][imem], markersize=7, mew = 0.5,
            marker=markers[imem], linewidth=0, color='r', label='Analysis Mem %s' % (markers[imem]) 
        )

    # Plot analysis increments
    for imem in range(ens_size):
        incre_x = all_data_dict[v1]['xa'][imem] - all_data_dict[v1]['xf'][imem]
        incre_y = all_data_dict[v2]['xa'][imem] - all_data_dict[v2]['xf'][imem]
        offset_x = incre_x
        offset_y = incre_y
        if imem==0:
            axs[1,icol].plot(
                    [all_data_dict[v1]['xf'][imem]+offset_x, all_data_dict[v1]['xa'][imem]-offset_x],
                    [all_data_dict[v2]['xf'][imem]+offset_y, all_data_dict[v2]['xa'][imem]-offset_y], 
                    color='r', linewidth=0.5,
                    label = 'Analysis Increment'
                )
        else:
            axs[1,icol].plot(
                    [all_data_dict[v1]['xf'][imem]+offset_x, all_data_dict[v1]['xa'][imem]-offset_x],
                    [all_data_dict[v2]['xf'][imem]+offset_y, all_data_dict[v2]['xa'][imem]-offset_y], 
                    color='r', linewidth=0.5
                )
    
    # Label plots
    axs[1,icol].set_xlabel( '%s (%s)' % (all_data_dict[v1]['name'], all_data_dict[v1]['unit']) )
    axs[1,icol].set_ylabel( '%s (%s)' % (all_data_dict[v2]['name'], all_data_dict[v2]['unit']) )
    axs[1,icol].set_title( 'b%d) Analysis %s & %s' % (icol+1, v2.upper(), v1.upper() ), loc='left')

    if v1 == 'qv':
        axs[1,icol].axvline( 2.5, color = 'k', zorder=1, linestyle="--", linewidth=0.5, label='Cloud-permitting QVAPOR')



# Throw on observation indicator
for i in range(2):
    for j in range(2):
        axs[i,j].axhline( obs, color='dodgerblue', linestyle='--', linewidth=1, label='Observed IR BT')



    


fig.subplots_adjust(wspace=0.4, hspace=0.5, left=0.075, right=0.99, bottom=0.25, top=0.93)


# Throw on some legends
handles, labels = axs[1,0].get_legend_handles_labels()
fig.legend( handles, labels, loc='lower center', ncols=4)


plt.savefig('figs/demo_dry_cloud_effect.pdf')

    
    




# # Plot QVAPOR VS IRBT FCST AND ANALYSIS MEMBERS
# # ---------------------------------------------

# # plot individual members
# for imem in range( len(irbt_xf_ens)):
#     axs[0,0].plot( irbt_xf_ens[imem], qvapor_xf_ens[imem], marker=markers[imem], linewidth=0, color='k', label='Forecast Mem %s' % (markers[imem]) )
#     axs[1,0].plot( irbt_xf_ens[imem], qvapor_xf_ens[imem], marker=markers[imem], linewidth=0, color='k', label='Forecast Mem %s' % (markers[imem]) )
#     axs[1,0].plot( irbt_xa_ens[imem], qvapor_xa_ens[imem], marker=markers[imem], linewidth=0, color='r', label='Analysis Mem %s' % (markers[imem]) )


# # Plot prior forecast PDF of QVAPOR and IRBT
# prior_cov = np.cov( irbt_xf_ens, qvapor_xf_ens )
# prior_mean = np.array( [np.mean(irbt_xf_ens), np.mean(qvapor_xf_ens)] )
# multivariate_normal 

# irvec = np.linspace( 205, 245, 101 )
# qvvec = np.linspace( 1.5,2.8, 102 )
# irmesh, qvmesh = np.meshgrid( irvec, qvvec)
# m


# # Plot labels


# axs[0].set_ylabel('6.2um BT (K)')
# axs[0].set_xlabel('QVAPOR (g/kg)')
# axs[0].set_title('a) Forecast BT & QVAPOR', loc='left')

# axs[1].set_ylabel('6.2um BT (K)')
# axs[1].set_xlabel('QHYDRO (g/kg)')
# axs[1].set_title('b) Forecast BT & QHYDROS', loc='left')



# # Plot QVAPOR VS IRBT 
# # --------------------

# # plot individual members
# for imem in range( len(irbt_xf_ens)):
#     axs[0].plot( qvapor_xf_ens[imem], irbt_xf_ens[imem], marker=markers[imem], linewidth=0, color='k', label='Forecast Mem %s' % (markers[imem]) )
# for imem in range( len(irbt_xf_ens)):
#     axs[0].plot( qvapor_xa_ens[imem], irbt_xa_ens[imem], marker=markers[imem], linewidth=0, color='r', label='Analysis Mem %s' % (markers[imem]) )

# # Plot increment lines
# for imem in range(len(irbt_xf_ens)):
#     incre_x = qvapor_xa_ens[imem] - qvapor_xf_ens[imem]
#     incre_y = irbt_xa_ens[imem] - irbt_xf_ens[imem]
#     offset_x = incre_x/np.sqrt( incre_x**2 + incre_y**2)
#     offset_y = incre_y/np.sqrt( incre_x**2 + incre_y**2)
#     if imem==0:
#         axs[0].plot(
#                 [qvapor_xf_ens[imem]+offset_x, qvapor_xa_ens[imem]-offset_x],
#                 [irbt_xf_ens[imem]+offset_y, irbt_xa_ens[imem]-offset_y], 
#                 color='r', linewidth=0.5,
#                 label = 'Analysis Increment'
#             )
#     else:
#         axs[0].plot(
#                 [qvapor_xf_ens[imem]+offset_x, qvapor_xa_ens[imem]-offset_x],
#                 [irbt_xf_ens[imem]+offset_y, irbt_xa_ens[imem]-offset_y], 
#                 color='r', linewidth=0.5
#             )






# # Plot QHYDRO VS IRBT 
# # --------------------

# # Plot ens members
# for imem in range( len(irbt_xf_ens)):
#     axs[1].plot( qhydro_xf_ens[imem], irbt_xf_ens[imem], marker=markers[imem], linewidth=0, color='k', label='Forecast Mem %s' % (markers[imem]) )
#     axs[1].plot( qhydro_xa_ens[imem], irbt_xa_ens[imem], marker=markers[imem], linewidth=0, color='r', label='Analysis Mem %s' % (markers[imem]) )

# # Plot increment lines
# for imem in range(len(irbt_xf_ens)):
#     incre_x = qhydro_xa_ens[imem] - qhydro_xf_ens[imem]
#     incre_y = irbt_xa_ens[imem] - irbt_xf_ens[imem]
#     offset_x = incre_x/np.sqrt( incre_x**2 + incre_y**2)
#     offset_y = incre_y/np.sqrt( incre_x**2 + incre_y**2)
#     axs[1].plot(
#         [ qhydro_xf_ens[imem]+offset_x, qhydro_xa_ens[imem]-offset_x],
#         [irbt_xf_ens[imem]+offset_y, irbt_xa_ens[imem]-offset_y], 
#         color='r', linewidth=0.5,
#         label = 'Analysis Increment'
#     )  





# # Throwing on some indicators of observations
# for i in range(2):
#     axs[i].axhline( obs, color='dodgerblue',zorder=0, linestyle='--', label='Observed BT' )

# # Throwing on saturation lines
# axs[0].axvline( 2.5, color = 'k', zorder=0, linestyle="--", linewidth=0.5, label='RH = 100%')

# fig.subplots_adjust(wspace=0.3, hspace=0.5, left=0.1, right=0.99, bottom=0.3, top=0.93)


# # Throw on some legends
# handles, labels = axs[0].get_legend_handles_labels()
# fig.legend( handles, labels, loc='lower center', ncols=4)

