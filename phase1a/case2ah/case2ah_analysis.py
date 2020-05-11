"""
PHASE 1A CASE 2AH 
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

sp = openmc.StatePoint('h5files/3pcm/fhr_p1a_c2ah_3pcm_statepoint.500.h5')
beta_b(sp,'c2ah')
fission_density_c(sp,'c2ah')
neutron_flux_d(sp,1.41823,0.00003,'c2ah')
neutron_flux_e(sp,1.41823,'c2ah')
neutron_spectrum_f(sp,'c2ah')