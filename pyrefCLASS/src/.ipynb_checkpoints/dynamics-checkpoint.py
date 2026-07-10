#from fncmb import ClassCosmo as Cosmo
import fos
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
plt.rc('text', usetex=True)

def make_timelist(initialz,z2eta,eta2z,Neta):
    ''' make conformal time list

    Parameters:
        Neta (int): stepsize (corresponding to Ntau in Nakatsuka code)
    Retrun: 
      
    '''

    initialeta = z2eta(initialz)
    finaleta = z2eta(0)
    etalist = np.linspace(initialeta,finaleta, Neta+1)
    #TODO：区切り方これでよいのか？最初のzが1e6になってないが
    etalist = np.delete(etalist,0)
    zlist = eta2z(etalist)

    return etalist, zlist

def load_dynamics(path_to_axionDynamics, etalist, zlist):
    ''' load axionDynamics.dat and transform it into DataFrame,
        adding dynamics of beta

    Parameters:
      path_to_axionDynamics (path)
      etalist
      zlist

    Retrun: 
      df (DataFrame): normalized phi, normalized beta from initial time, normalized beta from today, Note : mass(str) as df.index.meta 
    '''
    df = pd.read_table(path_to_axionDynamics,header=None)
    df = df.rename(columns={0:'phi'})
    df['beta(eta)'] = df['phi'].iloc[0] - df['phi']
    df['beta_0-beta(eta)'] = df['phi'].iloc[-1] - df['phi'] 
    df['eta'] = etalist
    df['z'] = zlist
    df = df[['eta','z','phi','beta(eta)','beta_0-beta(eta)']]

    filename = os.path.basename(path_to_axionDynamics)
    mass = fos.get_str(filename,'axionDynamics_m','_N')
    df.index.meta = {'mass':mass}
    
    return df

def plot_AxionDynamics(dfs):
    ''' plot normalized axion dynamics

    Parameters:
      dfs ([dataframe, dataframe, ...]) : list of dataframes of the return of load_dynamics func

    '''

    fig,ax = plt.subplots(figsize=(20,10),facecolor='white',edgecolor='white')

    for df in dfs:
        meta = df.index.meta
        mass = int(meta['mass'])
        zlist = df['z']
        philist = df['phi']
        ax.plot(zlist,philist,label=r'$\log_{10}m_{\phi}=-$'+f'{mass/100}'+' eV')
    
    ax.set_xscale('log')
    
    ax.set_xlabel('Redshift',fontsize='25')
    ax.set_ylabel(r'$\phi/\phi_{\mathrm{ini}}$',fontsize='35')

    ax.tick_params(axis='x', labelsize=35)
    ax.tick_params(axis='y', labelsize=35)

    # グラフにシャドウを追加
    reionization = patches.Rectangle((5, -50), 5, 60, linewidth=0, edgecolor=None, facecolor='lightblue', alpha=0.5)
    recombination = patches.Rectangle((868, -50), 376.7, 60, linewidth=0, edgecolor=None, facecolor='moccasin', alpha=0.9)

    plt.gca().add_patch(reionization)
    plt.gca().add_patch(recombination)
    
    plt.legend()
    ax.grid()

def plot_sin4beta_dynamics(dfs):
    ''' plot sin4beta dynamics

    Parameters:
      dfs ([dataframe, dataframe, ...]) : list of dataframes of the return of load_dynamics func

    '''

    fig,ax = plt.subplots(figsize=(20,10),facecolor='white',edgecolor='white')
    for df in dfs:
        meta = df.index.meta
        mass = int(meta['mass'])
        zlist = df['z']
        sin4betalist = np.sin(4*0.5*np.radians(0.34)*df['beta_0-beta(eta)'].values)
        ax.plot(zlist,sin4betalist,label=r'$\log_{10}m_{\phi}=-$'+f'{mass/100}'+' eV')

    ax.set_xscale('log')
    
    ax.set_xlabel('Redshift',fontsize='25')
    ax.set_ylabel(r'$\sin(4(\beta_0-\beta(z)))$',fontsize='35')

    ax.tick_params(axis='x', labelsize=35)
    ax.tick_params(axis='y', labelsize=35)

    # グラフにシャドウを追加
    reionization = patches.Rectangle((5, -50), 5, 60, linewidth=0, edgecolor=None, facecolor='lightblue', alpha=0.5)
    recombination = patches.Rectangle((868, -50), 376.7, 60, linewidth=0, edgecolor=None, facecolor='moccasin', alpha=0.9)

    plt.gca().add_patch(reionization)
    plt.gca().add_patch(recombination)

    plt.legend()
    ax.grid()

def plot_cos2beta_dynamics(dfs):
    ''' plot cos2beta dynamics
    Parameters:
      dfs ([dataframe, dataframe, ...]) : list of dataframes of the return of load_dynamics func

    '''

    fig,ax = plt.subplots(figsize=(20,10),facecolor='white',edgecolor='white')
    for df in dfs:
        meta = df.index.meta
        mass = int(meta['mass'])
        zlist = df['z']
        cos2betalist = np.cos(2*0.5*np.radians(0.34)*df['beta_0-beta(eta)'].values)
        ax.plot(zlist,cos2betalist,label=r'$\log_{10}m_{\phi}=-$'+f'{mass/100}'+' eV')

    ax.set_xscale('log')
    
    ax.set_xlabel('Redshift',fontsize='25')
    ax.set_ylabel(r'$\cos(2(\beta_0-\beta(z)))$',fontsize='35')

    ax.tick_params(axis='x', labelsize=35)
    ax.tick_params(axis='y', labelsize=35)

    # グラフにシャドウを追加
    reionization = patches.Rectangle((5, -50), 5, 60, linewidth=0, edgecolor=None, facecolor='lightblue', alpha=0.5)
    recombination = patches.Rectangle((868, -50), 376.7, 60, linewidth=0, edgecolor=None, facecolor='moccasin', alpha=0.9)

    plt.gca().add_patch(reionization)
    plt.gca().add_patch(recombination)

    plt.legend()
    ax.grid()
