"""
PHASE 1A CASE 1A 
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

sp = openmc.StatePoint('h5files/3pcm/fhr_p1a_c1a_3pcm_statepoint.500.h5')
beta_b(sp,'c1a')
fission_density_c(sp,'c1a')
neutron_flux_d(sp,1.40752,0.00003,'c1a')
neutron_flux_e(sp,1.40752,'c1a')
neutron_spectrum_f(sp,'c1a')