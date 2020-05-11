"""
PHASE 1A CASE 2AC
This python script builds analyzes the tallies from openmc statepoint file 
and uses functions in scripts/openmc_analysis.py to analyze and manipulate 
the data into what is required for the FHR benchmark. 
"""

###############################################################################
#                      Python Package Import
###############################################################################

import openmc
import sys
sys.path.insert(1, '../../scripts/')
from constants import *
from openmc_analysis import *

###############################################################################
#                                  Run
###############################################################################
case = '2ac'
keff = 1.43456
keff_unc = 0.00003

sp = openmc.StatePoint('h5files/3pcm/fhr_p1a_c2ac_3pcm_statepoint.500.h5')
beta_b(sp,case)
fission_density_c(sp,case)
neutron_flux_d(sp,keff,keff_unc,case)
neutron_flux_e(sp,keff,case)
neutron_spectrum_f(sp,case)