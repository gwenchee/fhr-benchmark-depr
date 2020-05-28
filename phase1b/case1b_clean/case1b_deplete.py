"""
This python script runs the depletion calculation for case1b of the FHR
benchmark.

"""

###############################################################################
#                      Python Package Import
###############################################################################
import openmc.deplete
from case1a_build_xml import *
import numpy as np

###############################################################################
#                                  Run
###############################################################################
chain_file = "../../data/chain_endfb71_pwr.xml"
chain = openmc.deplete.Chain.from_xml(chain_file)

operator = openmc.deplete.Operator(geom, settings, chain_file)



#time_steps = list(dep_time.copy() * 24 * 60 * 60)
hmop = operator.heavy_metal
print('hm = ' + str(hmop))

integrator = openmc.deplete.CELIIntegrator(operator=operator, timesteps=np.diff(bu), timestep_units='Mwd/kg', power=power_GW*1e9)
integrator.integrate()

