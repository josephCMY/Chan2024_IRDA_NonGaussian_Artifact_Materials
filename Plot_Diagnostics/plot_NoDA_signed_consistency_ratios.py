'''
    Script to examine the tails of the NoDA innovation distribution
'''

from scipy.stats import skewnorm, norm, moment
from matplotlib import use as mpl_use
mpl_use('agg')
import matplotlib.pyplot as plt
import numpy as np
from numba import njit
from gc import collect as gc_collect
from copy import deepcopy
from pylib_expt_specific_settings import *
from pylib_read_obs_priors import *

from datetime import datetime, timedelta

'''
    USER CONTROLS
'''
ens_size = 50
obs_err_sigma = 1.0

# Date controls
date_st = datetime.strptime('202301130000', '%Y%m%d%H%M')
date_ed = datetime.strptime('202301160000', '%Y%m%d%H%M')
time_interval = float( 60 )



'''
    Figure to show tail masses of IR BT
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

# Dictionary to hold onto CR values
cr_dict = {}
for vname in ['NoDA','Gauss', 'Std Norm']:
    cr_dict[vname] = []

# Dictionary to control plotting styles
plot_dict = {
    'NoDA': {'lw': 2, 'c':'r', 'zo': 10},
    'Gauss': {'lw': 2, 'c':'k', 'zo': 9},
    'Std Norm': {'lw': 4, 'c':'darkgray', 'zo': 8},
}



# Loop over all dates
for dd in range(nDates):

    date_nw = date_list[dd]

    # Load real data 
    fname = date_nw.strftime('DART_files/'+expt_name+'/obs_prior/obs_seq.final.%Y%m%d%H%M')
    f = open( fname, 'r' )
    obs_dict = read_obs_seq_prior_and_obsval( f, obs_kind=260, ens_size=ens_size )
    f.close()

    all_obs_prior_names = ['prior ensemble member     %2d' % (ee+1) \
                            for ee in range(ens_size)]
    all_obs_priors = np.array(
        [ obs_dict[all_obs_prior_names[ee]] for ee in range(ens_size) ] 
    )
    actual_obs = obs_dict[ 'observations']

    # Compute signed consistency ratio for individual obs & store in NoDA
    inno = actual_obs - np.mean( all_obs_priors, axis=0)
    sigma2 = (
        np.var( all_obs_priors, axis=0, ddof=1)
        + obs_err_sigma**2
    )
    cr_dict['NoDA'].append(
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
    cr_dict['Gauss'].append(
        np.sign(inno) * np.sqrt(
            (inno**2)/sigma2
        )
    )

    cr_dict['Std Norm'].append( np.random.normal( size = noda_prior_sigma.shape ) )

    print('Loaded ', date_nw )
# --- End of loop over dates




# Plot histogram of all masses
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




# For every 








