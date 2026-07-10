import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
plt.rc('text', usetex=True)
import utils


def plot_sin4beta_dynamics(df, makimakilist, xlim=None):
    ''' plot sin4beta dynamics for each makimaki parameter

    Parameters:
      dfs (dataframe) : dataframe of the return of dynamics.load_dynamics func
      makimakilist ([int, int,] or [float, float]) : list of makimakiparameter
      xlim ([min_redshift,max_redshift]) : set_xlim parameters
    '''
    fig, ax = plt.subplots(figsize=(20,10),facecolor='white',edgecolor='white')
    
    meta = df.index.meta
    mass = int(meta['mass'])
    zlist = df['z']

    for makimaki in makimakilist:
        sin4betalist = np.sin(4*0.5*np.radians(0.34+180*makimaki)*df['beta_0-beta(eta)'].values)
        ax.plot(zlist,sin4betalist,label=r'$m=$'+f'{makimaki}')

    ax.set_xscale('log')
    
    ax.set_xlabel('Redshift',fontsize='25')
    ax.set_ylabel(r'$\sin(4(\beta_0-\beta(z)))$',fontsize='35')

    if xlim is None:
        pass
    else:
      ax.set_xlim([xlim[0],xlim[1]])

    ax.tick_params(axis='x', labelsize=35)
    ax.tick_params(axis='y', labelsize=35)

    # グラフにシャドウを追加
    reionization = patches.Rectangle((5, -50), 5, 60, linewidth=0, edgecolor=None, facecolor='lightblue', alpha=0.5)
    recombination = patches.Rectangle((868, -50), 376.7, 60, linewidth=0, edgecolor=None, facecolor='moccasin', alpha=0.9)

    plt.gca().add_patch(reionization)
    plt.gca().add_patch(recombination)
    ax.text(5.2, -1, 'reionization', fontsize=25, color='black',rotation=90)
    ax.text(820, -1, 'recombination', fontsize=25, color='black',rotation=90)

    plt.title(r'$\log_{10}m_{\phi}=-$'+f'{mass/100}'+' eV',fontsize=25)
    plt.legend()
    ax.grid()


def plot_cos2beta_dynamics(df, makimakilist, xlim=None):
    ''' plot sin4beta dynamics for each makimaki parameter

    Parameters:
      dfs (dataframe) : dataframe of the return of dynamics.load_dynamics func
      makimakilist ([int, int,] or [float, float]) : list of makimakiparameter
      xlim ([min_redshift,max_redshift]) : set_xlim parameters
    '''
    fig, ax = plt.subplots(figsize=(20,10),facecolor='white',edgecolor='white')
    
    meta = df.index.meta
    mass = int(meta['mass'])
    zlist = df['z']

    for makimaki in makimakilist:
        cos2betalist = np.cos(2*0.5*np.radians(0.34+180*makimaki)*df['beta_0-beta(eta)'].values)
        ax.plot(zlist,cos2betalist,label=r'$m=$'+f'{makimaki}')

    ax.set_xscale('log')
    
    ax.set_xlabel('Redshift',fontsize='25')
    ax.set_ylabel(r'$\cos(2(\beta_0-\beta(z)))$',fontsize='35')

    if xlim is None:
        pass
    else:
      ax.set_xlim([xlim[0],xlim[1]])

    ax.tick_params(axis='x', labelsize=35)
    ax.tick_params(axis='y', labelsize=35)

    # グラフにシャドウを追加
    reionization = patches.Rectangle((5, -50), 5, 60, linewidth=0, edgecolor=None, facecolor='lightblue', alpha=0.5)
    recombination = patches.Rectangle((868, -50), 376.7, 60, linewidth=0, edgecolor=None, facecolor='moccasin', alpha=0.9)

    plt.gca().add_patch(reionization)
    plt.gca().add_patch(recombination)
    ax.text(5.2, -1, 'reionization', fontsize=25, color='black',rotation=90)
    ax.text(820, -1, 'recombination', fontsize=25, color='black',rotation=90)

    plt.title(r'$\log_{10}m_{\phi}=-$'+f'{mass/100}'+' eV',fontsize=25)
    plt.legend()
    ax.grid()


def get_Cl_by_maki(path_to_dir, files, correlation='EE'):
    ''' this func returns lists of power spectrm and makimaki parameter

    Parameters:
      path_to_dir (path) : path to directry of output files of birefCLASS
      files (list) : filenames which you want to use 
      correlation ('EE', 'EB', 'BB') : correlation which you want to plot
    Returns:
      ell (numpy) : muptipoleのnumpy配列
      PS_list (list) : a list of power spectrum
      legends_list (list) : makimaki parameters in a style of matplot legends 
    '''
    PS_dict={}
    for file in files:
        beta,_,_ = utils.get_basic_info(file)
        beta_deg = np.rad2deg(beta)
        maki = round(beta_deg/180)
        df = utils.load_output(path_to_dir + '/' + file)
        PS = df[correlation].values
        PS_dict[maki] = PS
    l = df['l'].values
    sorted_PS_dict = dict(sorted(PS_dict.items()))
    PS_list = list(sorted_PS_dict.values())
    maki_list = list(sorted_PS_dict.keys())
    legends_list = [f"${correlation}(m={value})$" for value in maki_list]
    return l,PS_list,legends_list


def get_Cl_by_delta(path_to_dir, files, correlation='EE'):
    ''' this func returns lists of power spectrm and delta parameter

    Parameters:
      path_to_dir (path) : path to directry of output files of birefCLASS
      files (list) : filenames which you want to use 
      correlation ('EE', 'EB', 'BB') : 
    
    Returns:
      ell (numpy) : muptipoleのnumpy配列
      PS_list (list) : a list of power spectrum
      legends_list (list) : delta parameters in a style of matplot legends 
    '''
    PS_dict={}
    for file in files:
        beta,_,_ = utils.get_basic_info(file)
        beta_deg = np.rad2deg(beta)
        delta = round(beta_deg - int(beta_deg/180)*180 -0.34, 2)
        df = utils.load_output(path_to_dir + '/' + file)
        PS = df[correlation].values
        PS_dict[delta] = PS
    l = df['l'].values
    sorted_PS_dict = dict(sorted(PS_dict.items()))
    PS_list = list(sorted_PS_dict.values())
    delta_list = list(sorted_PS_dict.keys())
    legends_list = [f"${correlation}$"+r" $(\delta$"+f"$={value})$" for value in delta_list]
    return l,PS_list,legends_list


