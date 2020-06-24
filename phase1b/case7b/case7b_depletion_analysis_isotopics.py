"""
This python script runs the isotopic depletion analysis for case7b of the FHR
benchmark.

"""
###############################################################################
#                      Python Package Import
###############################################################################

import openmc.deplete
import os
import pandas as pd
import sys
sys.path.insert(1, '../../scripts/')
from depletion_analysis import *
from openmc_analysis import *
from phase1b_constants import *

###############################################################################
#                                  Run
###############################################################################

case = 'p1b_c7b'
results = openmc.deplete.ResultsList.from_hdf5(
    "results/depletion_results_10000p_nobug.h5")
depletion_keff(results, 'long', case)
depletion_actinides(results, case, 'long')
depletion_fp(results, case, 'long')
depletion_extended(results, case, 'long')