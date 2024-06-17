'''
    Script to demonstrate that AOEI causes biases.
'''

from scipy.stats import skewnorm, norm, moment
from matplotlib import use as mpl_use
mpl_use('agg')
import matplotlib.pyplot as plt
import numpy as np
from numba import njit
from gc import collect as gc_collect

from pylib_read_obs_priors import *

from datetime import datetime, timedelta

'''
    USER CONTROLS
'''
ens_size = 50
prior_skewnorm_parameter = -10
num_trials = 3000

prior_mean, prior_var = skewnorm(prior_skewnorm_parameter).stats(moments='mv')
obs_err_sigma = 0.1

# Date controls
date_st = datetime.strptime('202301130000', '%Y%m%d%H%M')
date_ed = datetime.strptime('202301160000', '%Y%m%d%H%M')
time_interval = float( 60 )


'''
    GENERATE MASTER FIGURE
'''
fig, axs = plt.subplots( nrows=3, ncols=1, figsize=(6,6))


























'''
    Demonstrate that innovations involving ensemble means are skewed when the prior is skewed.
'''
prior_dist = skewnorm(prior_skewnorm_parameter)
# Draw ensemble and obs
fcst_ens = prior_dist.rvs( [num_trials, num_trials, ens_size] )
all_obs = prior_dist.rvs( [num_trials, num_trials] ) 
all_obs += np.random.normal(scale = obs_err_sigma, size=[num_trials, num_trials])



# ax.set_yscale('log')





# Plot out d = yo - <yf> distribution and mark AOEI, RAOEI and QC thresholds
# ---------------------------------------------------------------------------
ax = axs[1]
inno = all_obs - np.mean(fcst_ens, axis=2)
hist, edges = np.histogram( inno, bins=np.linspace(-4,4,101), density=True)
bin_ctrs = (edges[1:]+edges[:-1])/2

aoei_left_thresh = -np.sqrt(obs_err_sigma**2 + prior_var)
aoei_right_thresh = np.sqrt(obs_err_sigma**2 + prior_var)

raoei_left_thresh = -np.sqrt(obs_err_sigma**2 + prior_var)*3
raoei_right_thresh = np.sqrt(obs_err_sigma**2 + prior_var)*3

print( "Probability of CR < -3: ", np.mean( inno < -raoei_right_thresh) )
print( "Probability of CR > +3: ", np.mean( inno >  raoei_right_thresh) )

print( "Probability of CR < -1: ", np.mean( inno < -aoei_right_thresh) )
print( "Probability of CR > +1: ", np.mean( inno >  aoei_right_thresh) )

aoei_left_thresh_quantile = np.mean( inno < aoei_left_thresh )
aoei_right_thresh_quantile = 1-np.mean( inno < aoei_right_thresh )

cdf = np.cumsum(hist) / np.sum(hist)

ax.plot( bin_ctrs, hist, color='k', label = 'PDF of $\overline{d}$', zorder=2, linewidth=2)

ax.axvline( raoei_left_thresh, linestyle = '--', color='k', label='CR=3', zorder=0, linewidth=2)
ax.axvline( raoei_right_thresh, linestyle = '--', color='k', zorder=0, linewidth=2)

ax.axvline( aoei_left_thresh, linestyle = '--', color='r', label='AOEI thresh.', zorder=0, linewidth=2)
ax.axvline( aoei_right_thresh, linestyle = '--', color='r', zorder=0, linewidth=2)

print( 'mean inno of OmF', np.mean(inno),' +- ', np.std(inno, ddof=1)/num_trials)




# Plot out yo - yf_n distribution and mark AOEI thresholds
# --------------------------------------------------------
ax = axs[1]
inno = all_obs - ( fcst_ens[:,:,0]  )
hist, edges = np.histogram( inno, bins=np.linspace(-4,4,101), density=True)
bin_ctrs = (edges[1:]+edges[:-1])/2
aoei_left_thresh = -np.sqrt(obs_err_sigma**2 + 2*prior_var)
aoei_right_thresh = np.sqrt(obs_err_sigma**2 + 2*prior_var)
# ax.plot( bin_ctrs, hist, color='darkgray', label = 'PDF of $d_n \equiv y^o - y^f_n $', zorder=1, linewidth=0.5)
# ax.axvline( aoei_left_thresh, linestyle = '--', color='darkgray',  label='AOEI bounds for $d_n$', zorder=0, linewidth=0.5)
# ax.axvline( aoei_right_thresh, linestyle = '--', color='darkgray', zorder=0, linewidth=0.5)


# norm_dist = norm( loc=np.mean(inno), scale=np.std(inno, ddof=1) )
# ax.plot( bin_ctrs, norm_dist.pdf(bin_ctrs), linestyle='-', color='dodgerblue', zorder=0, label='Fitted Normal PDF')

aoei_left_thresh_quantile = np.mean( inno < aoei_left_thresh )
aoei_right_thresh_quantile = 1-np.mean( inno < aoei_right_thresh )
print( aoei_left_thresh_quantile, aoei_right_thresh_quantile)


# norm_dist = norm( loc=np.mean(inno), scale=np.std(inno, ddof=1) )
# ax.plot( bin_ctrs, norm_dist.pdf(bin_ctrs), linestyle='-', color='dodgerblue', zorder=0, label='Fitted Normal PDF')
ax.legend(loc='upper left')

# ax.set_yscale('log')

ax.axvline(0, color='lightgray', zorder=0)

print( aoei_left_thresh_quantile, aoei_right_thresh_quantile)

ax.set_xlim([-4,2])
ax.set_ylabel( 'PDF')
ax.set_ylim([1e-3, 1])
ax.set_yscale('log')
ax.set_xlabel( '$\overline{d}$ (thick black) or $d_n$ (thin gray)')
ax.set_title('b) Idealized Demo: PDF of $\overline{d}$ and $d_n$', loc='left')










'''
    Funcitons to sample incrmenets with different quality control schemes
'''

# Function to sample AOEI increment distribution for pair of prior and observation parameters. Assume 100 members
def sample_aoei_increments( prior_skew = 10., obs_err_std = 0.1, ens_size = 50, nTrials = 100 ):

    # Prior distribtuion
    prior_dist = skewnorm(prior_skew)
    prior_mean, prior_var = prior_dist.stats(moments='mv')

    # Draw obs and prior for all trials
    obs = prior_dist.rvs([nTrials, nTrials]) 
    obs += np.random.normal( scale=obs_err_std, size=nTrials)

    prior_ens = prior_dist.rvs([ens_size, nTrials, nTrials])

    
    
    # Generate prior mean and variance
    prior_mean = np.mean( prior_ens, axis=0 )
    prior_var = np.var( prior_ens, axis=0, ddof=1)

    # Compute increment under aoei
    inno = obs - prior_mean

    # Activate AOEI
    aoei_obs_var = np.zeros( [nTrials, nTrials] ) + obs_err_std**2
    flag_aoei = (inno**2 - prior_var) > obs_err_std**2
    aoei_obs_var[flag_aoei] = (inno**2 - prior_var)[flag_aoei]

    # print(inno.shape)
    incre = (prior_var / (aoei_obs_var + prior_var)) * inno
    # incre = (prior_var / (prior_var + obs_err_std**2) )*inno

    return incre







# Function to sample RAOEI increment distribution for pair of prior and observation parameters. Assume 100 members
def sample_raoei_increments( prior_skew = 10., obs_err_std = 0.1, ens_size = 50, nTrials = 100, threshold=3 ):

    # Prior distribtuion
    prior_dist = skewnorm(prior_skew)
    prior_mean, prior_var = prior_dist.stats(moments='mv')

    # Draw obs and prior for all trials
    obs = prior_dist.rvs([nTrials, nTrials]) 
    obs += np.random.normal( scale=obs_err_std, size=nTrials)

    prior_ens = prior_dist.rvs([ens_size, nTrials, nTrials])
    
    
    # Generate prior mean and variance
    prior_mean = np.mean( prior_ens, axis=0 )
    prior_var = np.var( prior_ens, axis=0, ddof=1)

    # Compute increment under aoei
    inno = obs - prior_mean

    # Activate RAOEI
    alpha2 = threshold **2
    raoei_obs_var = (inno**2 - alpha2*prior_var)/alpha2
    flag_no_infl = raoei_obs_var < obs_err_std**2
    raoei_obs_var[flag_no_infl] = obs_err_std**2

    # print(inno.shape)
    incre = (prior_var / (raoei_obs_var + prior_var)) * inno

    return incre










# Function to sample increments resulting from basic
def sample_gross_check_increments( prior_skew = 10., obs_err_std = 0.1, ens_size = 50, nTrials = 100, threshold=3 ):

    # Prior distribtuion
    prior_dist = skewnorm(prior_skew)
    prior_mean, prior_var = prior_dist.stats(moments='mv')

    # Draw obs and prior for all trials
    obs = prior_dist.rvs([nTrials, nTrials]) 
    obs += np.random.normal( scale=obs_err_std, size=nTrials)

    prior_ens = prior_dist.rvs([ens_size, nTrials, nTrials])
    
    # Generate prior mean and variance
    prior_mean = np.mean( prior_ens, axis=0 )
    prior_var = np.var( prior_ens, axis=0, ddof=1)

    # Compute increment 
    inno = obs - prior_mean

    # Activate QC
    sum_var = obs_err_std**2 + prior_var
    cr = np.sqrt(inno**2 / sum_var)
    
    # Otherwise, business as usual
    incre = (prior_var / (obs_err_std**2 + prior_var)) * inno

    flag_qc = cr > threshold
    incre[flag_qc] = 0

    return incre














# Estimate QC-related biases over a range of skewnorm parameters
skewparam_list = np.linspace(-10,10,21)
num_skewparams = len(skewparam_list)

# Init dictionary to hold outcomes
incre_bias_info = {}
for key in ['aoei','raoei','basic']:
    incre_bias_info[key] = {}
    for key2 in ['est','sig']:
        incre_bias_info[key][key2] = np.zeros(num_skewparams)

# Hold sampling functions
incre_bias_info['aoei']['func'] = sample_aoei_increments
incre_bias_info['raoei']['func'] = sample_raoei_increments
incre_bias_info['basic']['func'] = sample_gross_check_increments

# Actual experiments
for iskew in range(num_skewparams):

    for key in ['aoei','raoei','basic']:
        skewparam = skewparam_list[iskew]
        incre = incre_bias_info[key]['func']( prior_skew = skewparam, obs_err_std=obs_err_sigma, 
                                              ens_size = ens_size, nTrials=num_trials )
        incre_bias_info[key]['est'][iskew] = np.mean(incre)
        incre_bias_info[key]['sig'][iskew] = np.std(incre, ddof=1) / np.sqrt(np.prod(incre.shape)**2)

        gc_collect()

        print( 'Finished skewness parameter of ', skewparam, ' for qc ', key)










# Plot all increment biases and 2x std err range

# Colors
incre_bias_info['aoei']['color'] = 'r'
incre_bias_info['raoei']['color'] = 'dodgerblue'
incre_bias_info['basic']['color'] = 'k'

ax = axs[2]
for key in ['aoei','raoei','basic']:
    ax.plot( skewparam_list, incre_bias_info[key]['est'], 
             color = incre_bias_info[key]['color'], linewidth=2, 
             label = key.upper() )


ax.axhline(0, color='lightgray', zorder=0)
ax.axvline(0, color='lightgray', zorder=0)
ax.set_ylabel('Posterior Bias')
ax.set_xlabel('Skew Normal Shape Parameter (related to skewness)')
ax.set_title('c) Idealized Demo: QC-induced Posterior Bias', loc='left')
ax.legend(ncols=2)














'''
    Figure to show skewness of BTs
'''

# Expt to draw from
expt_name = 'NoDA'
ens_size=50


# Construct list of dates
date_list = []
date_nw = deepcopy( date_st ) 
while ( date_nw <= date_ed ):
    date_list.append( date_nw )
    date_nw += timedelta( minutes = time_interval )

# Determine number of dates
nDates = len( date_list )

skewness_list = [None for dd in range(nDates)]


# Loop over all dates
for dd in range(nDates):

    date_nw = date_list[dd]

    # Load real data 
    fname = date_nw.strftime('DART_files/'+expt_name+'/obs_prior/obs_seq.final.%Y%m%d%H%M')
    f = open( fname, 'r' )
    obs_dict = read_obs_seq_prior( f, obs_kind=260, ens_size=ens_size )
    f.close()

    all_obs_prior_names = ['prior ensemble member     %2d' % (ee+1) \
                            for ee in range(ens_size)]
    all_obs_priors = np.array(
        [ obs_dict[all_obs_prior_names[ee]] for ee in range(ens_size) ] 
    )

    # Estimate skewness
    third_moments = np.mean( 
        np.power(
            all_obs_priors - np.mean(all_obs_priors, axis=0), 
            3),
        axis=0 )
    third_moments *= ens_size**2 
    third_moments /= (ens_size-2)*(ens_size-1)

    skewness_list[dd] = third_moments / np.power( np.var( all_obs_priors, axis=0, ddof=1), 1.5)

    if dd == 0:
        print( np.percentile( skewness_list[dd],1), np.percentile( skewness_list[dd],99))
    
    print('Loaded ', date_nw )
# --- End of loop over dates


# Plot histogram of all third moments

all_prior_skewness = np.array(skewness_list)
skewness_bin_edges = np.linspace(-6,3,101)
skewness_hist, tmp = np.histogram( all_prior_skewness, bins = skewness_bin_edges, density=True)
skewness_cdf = np.cumsum( skewness_hist ) / np.sum(skewness_hist)

print('Quantile of 0 skewness:', np.mean( all_prior_skewness < 0 ) )

# Plot skewness statistics
ax = axs[0]
bin_ctrs = (skewness_bin_edges[1:]+skewness_bin_edges[:-1])/2
ax.plot( bin_ctrs, skewness_cdf, '-k', linewidth=2., label='CDF' )
ax.axvline(np.mean(all_prior_skewness), color='k', linestyle='-', label='Mean' )
ax.axvline(np.median(all_prior_skewness), color='k', linestyle='--', label='Median' )
ax.axvline(0, color='lightgrey', linestyle='-', zorder=0  )

ax.set_title('a) Distribution of 6.2 um BT Prior Skewness (NoDA Ensemble)', loc='left') 
ax.set_xlabel('Prior Skewness of 6.2 um BT (no units)')
ax.set_ylabel('CDF')
# ax.set_yscale('log')
ax.legend()


fig.subplots_adjust(top=0.95, bottom = 0.08, left=0.13, right=0.98, hspace=0.6)

plt.savefig('figs/demo_impacts_of_QC.pdf')

# --- End of loop over all dates




# For every 








