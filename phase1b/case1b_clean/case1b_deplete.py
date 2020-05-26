"""
This python script runs the depletion calculation for case1b of the FHR
benchmark.

"""

###############################################################################
#                      Python Package Import
###############################################################################
import openmc.deplete
from case1a_build_xml import *

###############################################################################
#                                  Run
###############################################################################

chain = openmc.deplete.Chain.from_xml("../../data/chain_endfb71_pwr.xml")

operator = openmc.deplete.Operator(geom, settings, "../../data/chain_endfb71_pwr.xml")



time_steps = list(dep_time.copy() * 24 * 60 * 60)
hmop = operator.heavy_metal
print('hm = ' + str(hmop))

integrator = openmc.deplete.CELIIntegrator(operator=operator, timesteps=time_steps, power_density=200)
integrator.integrate()

