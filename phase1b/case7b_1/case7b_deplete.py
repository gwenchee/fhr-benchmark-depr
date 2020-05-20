"""
This python script runs the depletion calculation for case7b of the FHR 
benchmark. 

"""

###############################################################################
#                      Python Package Import
###############################################################################
import openmc.deplete
from case4a_build_xml import * 

###############################################################################
#                                  Run
###############################################################################

chain = openmc.deplete.Chain.from_xml("../../data/chain_endfb71_pwr.xml")

operator = openmc.deplete.Operator(geom, settings, "../../data/chain_endfb71_pwr.xml")

time_steps = list(dep_time_7b.copy() * 24 * 60 * 60)

integrator = openmc.deplete.PredictorIntegrator(operator, time_steps, power_GW * 1e9)
integrator.integrate()
