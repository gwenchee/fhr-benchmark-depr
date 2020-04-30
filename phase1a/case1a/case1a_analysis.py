"""
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
#                      Run
###############################################################################

sp = openmc.StatePoint('statepoint.10.h5')
print(beta_b(beta,sp))
fission_density_c(sp,'c1a')
neutron_flux_d(sp,1.47061,0.003776,'c1a')
neutron_flux_e(sp,1.47061,'c1a')
neutron_spectrum_f(sp,'c1a')