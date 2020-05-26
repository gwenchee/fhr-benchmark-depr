""" Functions for analyzing openmc tally results 

This scripts contains functions for analyzing the tallies 
from the openmc statepoint file and manipulate the data 
into what is required for the FHR benchmark. 

"""

import numpy as np
import pylab as pl
import matplotlib.colorbar as cbar
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.insert(1, '../../scripts/')
from phase1a_constants import *

###############################################################################
#                           Criticality Functions
###############################################################################


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


def beta_b(sp,case): 
    """Returns Beta-effective and its uncertainty

    Parameters
    ----------
    sp: openmc.statepoint.StatePoint
        this statepoint.h5 file is created by running the openmc 
        executable on the xml files generated by the build_xml.py 
        files and when tallies_on toggle is on True. 
    case: str 
        case number for naming files 
    Returns
    -------
    This function generates a csv file with beta-effective 
    and its uncertainty.
    """

    name = 'analysis_output/p1a_'+case+'_b'
    mesh_tally_b = sp.get_tally(name='mesh tally b')
    df_b = mesh_tally_b.get_pandas_dataframe()
    beta = df_b['mean'][0]/df_b['mean'][1]
    beta_err = beta * np.sqrt((df_b['std. dev.'][0]/df_b['mean'][0])**2+(df_b['std. dev.'][1]/df_b['mean'][1])**2)
    df_bb = pd.DataFrame()
    df_bb['beta'] = [beta]
    df_bb['beta err'] = [beta_err]
    df_bb.to_csv(name+'.csv')
    return 


def reactivity_coefficient_b(keff_og,keff_og_unc,keff_new,keff_new_unc,temp_change):
    """Generates the reactivity coefficient and its uncertainty

    Parameters
    ----------
    keff_og: float 
        original keff 
    keff_og_unc: float 
        original keff's uncertainty 
    keff_new: float 
        keff after temperature change
    keff_og_new: float 
        keff's uncertainty after temperature change
    temp_change: float 
        temperature change (be sure to include +/- sign)

    Returns
    -------
    coeff: float 
        reactivity coefficient
    coeff_unc: float 
        reactivity coefficient uncertainty
    """

    coeff = (keff_new*1e5 - keff_og*1e5) / temp_change
    coeff_unc = np.sqrt((keff_og_unc*1e5)**2 + (keff_new_unc*1e5)**2) / temp_change
    
    return coeff, coeff_unc


def fission_density_c(sp,case): 
    """Generates a csv and png file with results of fission source
    distribution by 1/5 stripes for phase 1a-c of the benchmark
    in the analysis_output folder. 

    Parameters
    ----------
    sp: openmc.statepoint.StatePoint
        this statepoint.h5 file is created by running the openmc 
        executable on the xml files generated by the build_xml.py 
        files and when tallies_on toggle is on True. 
    case: str 
        case number for naming files 
    Returns
    -------
    This function generates a csv file with fission density results 
    and visualization of fission source distribution by 1/5 stripe
    for phase 1a-c of the benchmark. 
    """

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
    ave_sd = 1/(len(df['mean'])) * np.sqrt(np.sum(df['std. dev.']**2))
    df['FD std dev'] = df['Fission Density']* np.sqrt((df['std. dev.']/df['mean'])**2 + (ave_sd/ave)**2)
    df['Relative unc.'] =  df['FD std dev']/df['Fission Density']
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
    """Generates a csv file with results of neutron flux
    averaged over the whole model, tabulated in 3 coarse energy groups 
    (upper energy boundaries 3 eV for thermal group and 0.1 MeV for 
    intermediate group)

    Parameters
    ----------
    sp: openmc.statepoint.StatePoint
        this statepoint.h5 file is created by running the openmc 
        executable on the xml files generated by the build_xml.py 
        files and when tallies_on toggle is on True. 
    k: float
        k-effective 
    kerr: float 
        k-effective uncertainty 
    case: str 
        case number for naming files 
    Returns
    -------
    This function generates a csv file with neutron flux results 
    """

    name = 'analysis_output/p1a_'+case+'_d'
    mesh_tally_d = sp.get_tally(name='mesh tally d')
    df_d = mesh_tally_d.get_pandas_dataframe()
    df_dd = pd.DataFrame(index=['E3','E2','E1'])
    df_dd['flux'], df_dd['flux_err'] = flux_conv(df_d,200,k,kerr)
    df_dd['relative_err_p'] = df_dd['flux_err'] / df_dd['flux'] 
    df_dd = df_dd.reindex(['E1','E2','E3'])
    df_dd.to_csv(name+'.csv')
    return 


def neutron_flux_e(sp,k,case): 
    """Generates a csv file and png files with results of neutron flux
    at 10000 points in model, tabulated in 3 coarse energy groups 
    (upper energy boundaries 3 eV for thermal group and 0.1 MeV for 
    intermediate group)

    Parameters
    ----------
    sp: openmc.statepoint.StatePoint
        this statepoint.h5 file is created by running the openmc 
        executable on the xml files generated by the build_xml.py 
        files and when tallies_on toggle is on True. 
    k: float
        k-effective 
    case: str 
        case number for naming files 
    Returns
    -------
    This function generates a csv file with neutron flux results 
    at 10000 points in the model for 3 energy groups, and 3 png 
    files visualizing the neutron flux distribution for 3 
    energy groups.  
    """

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
    plt.imshow(flux_conv['eg1']/np.mean(flux_conv['eg1']), interpolation='none', origin='lower',cmap='viridis')
    plt.colorbar()
    plt.title('Case '+case+' Energy Group 1 Flux Distribution')
    plt.savefig(name+'_eg1')
    np.savetxt(name+"_eg1.csv", np.flip(flux_conv['eg1'],0), delimiter=",")

    plt.figure()
    plt.imshow(flux_conv['eg2']/np.mean(flux_conv['eg2']), interpolation='none', origin='lower',cmap='viridis')
    plt.colorbar()
    plt.title('Case '+case+ ' Energy Group 2 Flux Distribution')
    plt.savefig(name+'_eg2')
    np.savetxt(name+"_eg2.csv", np.flip(flux_conv['eg2'],0), delimiter=",")

    plt.figure()
    plt.imshow(flux_conv['eg3']/np.mean(flux_conv['eg3']), interpolation='none', origin='lower',cmap='viridis')
    plt.colorbar()
    plt.title('Case '+case+ ' Energy Group 3 Flux Distribution')
    plt.savefig(name+'_eg3')
    np.savetxt(name+"_eg3.csv", np.flip(flux_conv['eg3'],0), delimiter=",")
    return 


def neutron_spectrum_f(sp,case,k,kerr): 
    """Generates a csv file and png file with results of neutron
    spectrum averaged over the fuel assembly. 

    Parameters
    ----------
    sp: openmc.statepoint.StatePoint
        this statepoint.h5 file is created by running the openmc 
        executable on the xml files generated by the build_xml.py 
        files and when tallies_on toggle is on True. 
    case: str 
        case number for naming files 
    Returns
    -------
    This function generates a csv and png file with results of neutron
    spectrum averaged over the fuel assembly. 
    """

    name = 'analysis_output/p1a_'+case+'_f'
    mesh_tally_f = sp.get_tally(name='mesh tally f')
    df_f = mesh_tally_f.get_pandas_dataframe()
    index_list = []
    for x in range(252):
        index_list += ['E'+str(x+1)]
        
    df_ff = pd.DataFrame(index=index_list)

    df_ff['flux'], df_ff['flux_err'] = flux_conv(df_f,200,k,kerr)
    fluxvals = np.append(np.array(df_ff['flux']),np.array(df_ff['flux'])[0])
    plt.figure()
    plt.semilogx(np.array(engs)/1e6,fluxvals,drawstyle='steps')
    plt.title('Case '+case+ ' Average neutron spectrum')
    plt.xlabel('energy [MeV]')
    plt.ylabel('flux [$n/cm^2s$]')
    plt.ylim(0, 1.3e13)
    plt.savefig(name)
    df_ff_T = df_ff.T
    df_ff_T.to_csv(name+'.csv')
    return 
