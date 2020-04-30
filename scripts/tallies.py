import openmc 
import openmc.mgxs as mgxs
from constants import *

def tallies_generation(root): 
    tallies_file = openmc.Tallies()
    energy_groups = openmc.mgxs.EnergyGroups(group_edges=[1e-6, 20.0e6])
    delayed_groups = list(range(1,7))
    # phase1a-b
    beta = mgxs.Beta(domain=root, energy_groups=energy_groups, delayed_groups=delayed_groups)
    beta.nuclides = ['U235']
    tallies_file += beta.tallies.values()
    # phase1a-c 
    mesh_no = 0
    for t in range(6):
        mesh_no += 1
        for x in range(2):
            x_trans = t*T['A1']['P']['x']
            y_trans = t*T['A1']['P']['y']
            if x == 1:
                mesh_no += 1
                x_trans += T['A1']['F']['x']
                y_trans += T['A1']['F']['y']
            mesh_c = openmc.RegularMesh(mesh_id=mesh_no)
            mesh_c.dimension = [1,5]
            mesh_c.lower_left = [V['A1']['F']['L']['x']+x_trans, 
                            V['A1']['F']['B']['y']+y_trans]
            mesh_c.upper_right = [V['A1']['F']['R']['x']+x_trans, 
                                V['A1']['F']['T']['y']+y_trans]
            mesh_filter_c = openmc.MeshFilter(mesh_c)
            tally_c = openmc.Tally(name='mesh tally c'+str(mesh_no))
            tally_c.filters = [mesh_filter_c]
            tally_c.scores = ['fission']
            tallies_file.append(tally_c)    
    # phase 1a-d
    energy_filter_d = openmc.EnergyFilter([1e-5, 3, 1.0e5,20.0e6])
    mesh_d = openmc.RegularMesh(mesh_id=13)
    mesh_d.dimension = [1,1]
    L = 27.02
    mesh_d.lower_left = [-L,-L]
    mesh_d.upper_right = [L,L]
    mesh_filter_d = openmc.MeshFilter(mesh_d)
    tally_d = openmc.Tally(name='mesh tally d')
    tally_d.filters = [mesh_filter_d, energy_filter_d]
    tally_d.scores = ['flux','nu-fission','fission']
    tallies_file.append(tally_d)
    # phase 1a-e
    energy_filter_e = openmc.EnergyFilter([1e-5, 3, 0.1e6,20.0e6])
    mesh_e = openmc.RegularMesh(mesh_id=14)
    mesh_e.dimension = [100,100]
    L = 27.02
    mesh_e.lower_left = [-L,-L]
    mesh_e.upper_right = [L,L]
    mesh_filter_e = openmc.MeshFilter(mesh_e)
    tally_e = openmc.Tally(name='mesh tally e')
    tally_e.filters = [mesh_filter_e, energy_filter_e]
    tally_e.scores = ['flux','nu-fission','fission']
    tallies_file.append(tally_e)
    # phase 1a-f 
    engs = [1.00E-11,1.00E-10,5.00E-10,7.50E-10,1.00E-09,1.20E-09,
            1.50E-09,2.00E-09,2.50E-09,3.00E-09,4.00E-09,5.00E-09,
            7.50E-09,1.00E-08,2.53E-08,3.00E-08,4.00E-08,5.00E-08,
            6.00E-08,7.00E-08,8.00E-08,9.00E-08,1.00E-07,1.25E-07,
            1.50E-07,1.75E-07,2.00E-07,2.25E-07,2.50E-07,2.75E-07,
            3.00E-07,3.25E-07,3.50E-07,3.75E-07,4.00E-07,4.50E-07,
            5.00E-07,5.50E-07,6.00E-07,6.25E-07,6.50E-07,7.00E-07,
            7.50E-07,8.00E-07,8.50E-07,9.00E-07,9.25E-07,9.50E-07,
            9.75E-07,1.00E-06,1.01E-06,1.02E-06,1.03E-06,1.04E-06,
            1.05E-06,1.06E-06,1.07E-06,1.08E-06,1.09E-06,1.10E-06,
            1.11E-06,1.12E-06,1.13E-06,1.14E-06,1.15E-06,1.18E-06,
            1.20E-06,1.23E-06,1.25E-06,1.30E-06,1.35E-06,1.40E-06,
            1.45E-06,1.50E-06,1.59E-06,1.68E-06,1.77E-06,1.86E-06,
            1.94E-06,2.00E-06,2.12E-06,2.21E-06,2.30E-06,2.38E-06,
            2.47E-06,2.57E-06,2.67E-06,2.77E-06,2.87E-06,2.97E-06,
            3.00E-06,3.10E-06,3.20E-06,3.50E-06,3.73E-06,4.10E-06,
            4.70E-06,5.00E-06,5.40E-06,6.00E-06,6.25E-06,6.50E-06,
            6.75E-06,6.88E-06,7.00E-06,7.15E-06,8.10E-06,9.10E-06,
            1.00E-05,1.15E-05,1.19E-05,1.29E-05,1.44E-05,1.60E-05,
            1.70E-05,1.85E-05,1.94E-05,2.00E-05,2.05E-05,2.12E-05,
            2.18E-05,2.25E-05,2.50E-05,2.75E-05,3.00E-05,3.13E-05,
            3.18E-05,3.33E-05,3.38E-05,3.50E-05,3.55E-05,3.60E-05,
            3.70E-05,3.71E-05,3.73E-05,3.76E-05,3.80E-05,3.91E-05,
            3.96E-05,4.10E-05,4.24E-05,4.40E-05,4.52E-05,4.83E-05,
            5.06E-05,5.34E-05,5.80E-05,6.10E-05,6.30E-05,6.50E-05,
            6.75E-05,7.20E-05,7.60E-05,8.00E-05,8.17E-05,9.00E-05,
            9.70E-05,1.01E-04,1.05E-04,1.08E-04,1.13E-04,1.16E-04,
            1.18E-04,1.19E-04,1.22E-04,1.43E-04,1.70E-04,1.80E-04,
            1.88E-04,1.89E-04,1.92E-04,1.93E-04,2.02E-04,2.07E-04,
            2.10E-04,2.20E-04,2.40E-04,2.85E-04,3.05E-04,5.50E-04,
            6.70E-04,6.83E-04,9.50E-04,1.15E-03,1.50E-03,1.55E-03,
            1.80E-03,2.20E-03,2.25E-03,2.50E-03,3.00E-03,3.74E-03,
            3.90E-03,5.70E-03,8.03E-03,9.50E-03,1.30E-02,1.70E-02,
            2.00E-02,3.00E-02,4.50E-02,5.00E-02,5.20E-02,6.00E-02,
            7.30E-02,7.50E-02,8.20E-02,8.50E-02,1.00E-01,1.28E-01,
            1.49E-01,2.00E-01,2.70E-01,3.30E-01,4.00E-01,4.20E-01,
            4.40E-01,4.70E-01,4.92E-01,5.50E-01,5.73E-01,6.00E-01,
            6.70E-01,6.79E-01,7.50E-01,8.20E-01,8.61E-01,8.75E-01,
            9.00E-01,9.20E-01,1.01E+00,1.10E+00,1.20E+00,1.25E+00,
            1.32E+00,1.36E+00,1.40E+00,1.50E+00,1.85E+00,2.35E+00,
            2.48E+00,3.00E+00,4.30E+00,4.80E+00,6.43E+00,8.19E+00,
            1.00E+01,1.28E+01,1.38E+01,1.46E+01,1.57E+01,1.73E+01,
            2.00E+01]
    engs = [x*1e6 for x in engs]
    energy_filter_f = openmc.EnergyFilter(engs)
    mesh_f = openmc.RegularMesh(mesh_id=15)
    mesh_f.dimension = [1,1]
    L = 27.02
    mesh_f.lower_left = [-L,-L]
    mesh_f.upper_right = [L,L]
    mesh_filter_f = openmc.MeshFilter(mesh_f)
    tally_f = openmc.Tally(name='mesh tally f')
    tally_f.filters = [mesh_filter_f, energy_filter_f]
    tally_f.scores = ['flux','nu-fission','fission']
    tallies_file.append(tally_f)
    
    tallies_file.export_to_xml()
    return beta, engs