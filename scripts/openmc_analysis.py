""" Functions for analyzing openmc tally results 

This scripts contains functions for analyzing the tallies 
from the openmc statepoint file and manipulate the data 
into what is required for the FHR benchmark. 

"""

#from pyne import data
import numpy as np
import pylab as pl
import matplotlib.colorbar as cbar
import pandas as pd
import matplotlib.pyplot as plt
from case1a_build_xml import * 
import sys
sys.path.insert(1, '../../scripts/')
from constants import *


###############################################################################
#                                  Functions
###############################################################################


def reactor_power(sp_power, particles):
    """Returns the FHR's power 

    Parameters
    ----------
    sp_power: float 
        specific power of reactor [W/gU]
    particles: int 
        number of triso particles in each fuel stripes's y-direction

    Returns
    -------
    t_u: float 
        mass of Uranium [metric tonnes]
    power: float 
        power of reactor [W]
    """

    wtp_u = ((0.00227325*data.atomic_mass('U235'))+ \
            (0.02269476*data.atomic_mass('U238')))/ \
            (0.00227325*data.atomic_mass('U235')+ \
             0.02269476*data.atomic_mass('U238')+ \
             0.03561871*data.atomic_mass('O16') + 
             0.00979714*data.atomic_mass('C0'))*100 #%
    V = (101 * 210 * particles * 36 * 4/3 * np.pi * (2135e-5)**3) #cm3
    density = 11 #g/cc
    grams_u = V * density * (wtp_u) / 100 #gU
    t_u = grams_u * 1e-6 #t
    power = grams_u * sp_power # W

    return t_u, power 


def flux_conv(df,sp_power,k,kerr):
    """Converts flux unit from [n*cm/src] to [n/cm^2*s]

    Parameters
    ----------
    df: pandas dataframe with flux, nu-fission, and fission values
    sp_power: float 
        specific power of reactor [W/gU]
    k: keff [n/src]
    kerr: keff uncertainty [n/src]

    Returns
    -------
    flux: np.array(float)
        array of flux in [n/cm2*s] for each energy group
    flux_err: np.array(float)
        array of flux uncertainties in [n/cm2*s] for each energy group

    """
    P = 245486.6796001383 # W 
    Q = 200*1.6022e-13 # J/fission
    nu_fission = np.array(df[df['score'].str.match('nu-fission')]['mean']) # n/src
    fission = np.array(df[df['score'].str.match('fission')]['mean']) # fission/src
    og_flux = np.array(df[df['score'].str.match('flux')]['mean']) # n*cm/src
    nu_fission_err = np.array(df[df['score'].str.match('nu-fission')]['std. dev.'])
    fission_err = np.array(df[df['score'].str.match('fission')]['std. dev.'])
    og_flux_err = np.array(df[df['score'].str.match('flux')]['std. dev.'])
    nu = nu_fission/fission # n/fission
    N = P*nu/(Q*k) # src/s
    V = 3 * np.sqrt(3) / 2 * H_side ** 2 * z_thickness * T_pitch # cm3
    flux = 1/ V * N * og_flux # n/(cm2*s)
    flux[np.isnan(flux)] = 0
    flux_err = (np.sqrt((nu_fission_err/nu_fission)**2+\
               (fission_err/fission)**2+(og_flux_err/og_flux)**2 + \
                (kerr/k)**2)) * flux
    flux_err[np.isnan(flux_err)] = 0

    return flux, flux_err 


def beta_b(beta,sp): 
    """Returns Beta-effective

    Parameters
    ----------
    beta: openmc.mgxs.mdgxs.Beta
    sp: openmc.statepoint.StatePoint

    Returns
    -------
    beff: float 
        mass of Uranium [metric tonnes]
    beff_err: float 
        power of reactor [W]
    """

    beta.load_from_statepoint(sp)
    beta.build_hdf5_store(filename='mgxs', append=True)
    df = beta.get_pandas_dataframe()
    beff = df['mean'].sum()
    beff_err = np.sqrt((df['std. dev.']**2).sum())
    return beff, beff_err


def fission_density_c(sp,case): 
    name = 'analysis_output/p1a_'+case+'_c'
    region = ['1','2','3','4','5']
    fission_rates = []
    num = 1
    for x in range(1,13): 
        mesh_tally = sp.get_tally(name='mesh tally c'+str(x))
        a = mesh_tally.get_pandas_dataframe()
        a = a.drop(columns=['mesh '+str(x),'nuclide','score'])
        if (x%2) == 0:
            num = int(x/2)
            stripe = [str(num)+'B']*5
        else:
            num = int((x+1)/2)
            stripe = [str(num)+'T']*5
        a['Region'] = region
        a['Stripe'] = stripe
        if x == 1: 
            df = a
        else:
            df = df.append(a)
        b = a['mean'].to_numpy()
        fission_rates.append(b)
    df = df.set_index(['Stripe','Region'])
    ave = df['mean'].mean()
    df['Fission Density'] = df['mean']/ave
    ave_sd = np.sqrt(((df['mean'] - ave)**2).mean())
    df['FD std dev'] = np.sqrt((df['std. dev.']/df['mean'])**2 + (ave_sd/ave)**2)
    df['Relative unc.'] =  df['FD std dev']/df['Fission Density']*100
    df.to_csv(name+'.csv')

    xs = []
    ys = []
    ws = np.array([F_len/5]*60)
    hs = np.array([F_width]*60)
    fission_rates /= np.mean(fission_rates)
    vs = fission_rates.flatten()

    for p in range(6): 
        for f in range(2): 
            x_trans = p*T['A1']['P']['x']
            y_trans = p*T['A1']['P']['y']
            if f == 1:
                x_trans += T['A1']['F']['x']
                y_trans += T['A1']['F']['y']
            for s in range(5):
                if s > 0:
                    x_trans += F_len/5
                xs.append(V['A1']['F']['L']['x']+x_trans)
                ys.append(V['A1']['F']['B']['y']+y_trans)

    normal = pl.Normalize(vs.min(), vs.max())
    colors = pl.cm.YlOrRd(normal(vs))

    ax = pl.subplot(111)
    for x,y,w,h,c in zip(xs,ys,ws,hs,colors):
        rect = pl.Rectangle((x,y),w,h,color=c)
        ax.add_patch(rect)

    cax, _ = cbar.make_axes(ax) 
    cb2 = cbar.ColorbarBase(cax, cmap=pl.cm.YlOrRd,norm=normal) 

    ax.set_xlim(-25,12)
    ax.set_ylim(-25,0)
    ax.set_xlabel('x [cm]')
    ax.set_ylabel('y [cm]')
    ax.set_title('Case ' + case + ': Normalized Fission Source')
    pl.savefig(name,bbox_inches='tight')
    return 


def neutron_flux_d(sp,k,kerr,case): 
    name = 'analysis_output/p1a_'+case+'_d'
    mesh_tally_d = sp.get_tally(name='mesh tally d')
    df_d = mesh_tally_d.get_pandas_dataframe()
    df_dd = pd.DataFrame(index=['E3','E2','E1'])
    df_dd['flux'], df_dd['flux_err'] = flux_conv(df_d,200,k,kerr)
    df_dd['relative_err_p'] = df_dd['flux_err'] / df_dd['flux'] * 100
    df_dd = df_dd.reindex(['E1','E2','E3'])
    df_dd.to_csv(name+'.csv')
    return 


def neutron_flux_e(sp,k,case): 
    name = 'analysis_output/p1a_'+case+'_e'
    mesh_tally_e = sp.get_tally(name='mesh tally e')
    flux = mesh_tally_e.get_slice(scores=['flux'])
    nu_fission = mesh_tally_e.get_slice(scores=['nu-fission'])
    fission = mesh_tally_e.get_slice(scores=['fission'])
    flux_conv = {}
    eg_names = ['eg3','eg2','eg1']
    egs = [(1e-5,3),(3,0.1e6),(0.1e6,20e6)]
    P = 245486.6796001383
    Q = 200*1.6022e-13
    V = 3 * np.sqrt(3) / 2 * H_side ** 2 * z_thickness * T_pitch/(100*100)
    for x in range(3): 
        flux_eg = flux.get_slice(filters=[openmc.EnergyFilter], filter_bins=[(egs[x],)])
        nu_fiss_eg = nu_fission.get_slice(filters=[openmc.EnergyFilter], filter_bins=[(egs[x],)])
        fiss_eg = fission.get_slice(filters=[openmc.EnergyFilter], filter_bins=[(egs[x],)])
        nu = nu_fiss_eg.mean / fiss_eg.mean
        nu = np.nanmean(nu)
        N = P*nu/(Q*k)
        flux_conv[eg_names[x]] = flux_eg.mean * 1/ V * N 
        flux_conv[eg_names[x]].shape = (100,100)
        flux_conv[eg_names[x]][np.isnan(flux_conv[eg_names[x]])] = 0

    plt.figure()
    plt.imshow(flux_conv['eg1'], interpolation='none', origin='lower')
    plt.colorbar()
    plt.title('Energy Group 1 Flux Distribution')
    plt.savefig(name+'_eg1')
    np.savetxt(name+"_eg1.csv", np.flip(flux_conv['eg1'],0), delimiter=",")

    plt.figure()
    plt.imshow(flux_conv['eg2'], interpolation='none', origin='lower')
    plt.colorbar()
    plt.title('Energy Group 2 Flux Distribution')
    plt.savefig(name+'_eg2')
    np.savetxt(name+"_eg2.csv", np.flip(flux_conv['eg2'],0), delimiter=",")

    plt.figure()
    plt.imshow(flux_conv['eg3'], interpolation='none', origin='lower')
    plt.colorbar()
    plt.title('Energy Group 3 Flux Distribution')
    plt.savefig(name+'_eg3')
    np.savetxt(name+"_eg3.csv", np.flip(flux_conv['eg3'],0), delimiter=",")
    return 


def neutron_spectrum_f(sp,case): 
    name = 'analysis_output/p1a_'+case+'_f'
    mesh_tally_f = sp.get_tally(name='mesh tally f')
    df_f = mesh_tally_f.get_pandas_dataframe()
    index_list = []
    for x in range(252):
        index_list += ['E'+str(x+1)]
        
    df_ff = pd.DataFrame(index=index_list)

    df_ff['flux'], df_ff['flux_err'] = flux_conv(df_f,200,1.42807,0.003776)
    fluxvals = np.append(np.array(df_ff['flux']),np.array(df_ff['flux'])[0])
    plt.figure()
    plt.semilogx(np.array(engs)/1e6,fluxvals,drawstyle='steps')
    plt.title('Average neutron spectrum')
    plt.xlabel('energy [MeV]')
    plt.ylabel('flux [$n/cm^2s$]')
    plt.savefig(name)
    df_ff_T = df_ff.T
    df_ff_T.to_csv(name+'.csv')
    return 