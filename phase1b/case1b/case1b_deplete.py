"""
This python script runs the depletion calculation for case1b of the FHR 
benchmark. 

"""

###############################################################################
#                      Python Package Import
###############################################################################

import openmc.deplete
import sys
sys.path.insert(1, '../../phase1a/case1a/')
from case1a_build_xml import * 

###############################################################################
#                                  Run
###############################################################################

chain = openmc.deplete.Chain.from_xml("../../data/chain_endfb71_pwr.xml")

operator = openmc.deplete.Operator(geom, settings, "../../data/chain_endfb71_pwr.xml")

time_steps = list(dep_time.copy() * 24 * 60 * 60)[:2]

integrator = openmc.deplete.PredictorIntegrator(operator, time_steps, power_GW * 1e9)
integrator.integrate()