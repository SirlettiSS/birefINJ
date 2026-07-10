import fos
import os
import numpy as np
import pandas as pd

def get_basic_info(filename):
    ''' extract beta(rad), mass(eg:2800), neta from output file name.

    Parameters:
      filename (str) : file name of birefCLASS output

    returns:
      beta (float) : birefringence angle in radian
      mass (int) : mass of ALP (10^{-31.2}eV => 3120)
      neta (int) : stepsize of etalist
    '''
    beta_str = fos.get_str(filename,'alpha0','_m')
    mass_str = fos.get_str(filename,'_m','_N')
    if '+' in filename:
        neta_str = fos.get_str(filename, '_N', '+cl')
    elif '-' in filename:
        neta_str = fos.get_str(filename, '_N', '-cl')
    else:
        print('Error')
    return int(beta_str)*1e-4,int(mass_str),int(neta_str)

def cp_output(path_to_output_birefCLASS, path_to_dir_save, time_crit):
    ''' ある基準時以降に作られたbirefCLASSのoutputのうち、cl_lensedのみをコピーします

    Paramters:
      path_to_dir_output_birefCLASS (path) : birefCLASS の output ディレクトリ
      path_to_dir_save (path) : コピー先のディレクトリ
      datetime_crit (datetime) : 基準時刻
    '''  
    if not os.path.exists(path_to_dir_save):
        os.mkdir(path_to_dir_save)
    fos.cp_by_timecriteria(path_to_output_birefCLASS,path_to_dir_save,time_crit)
    fos.select(path_to_dir_save, ['cl_lensed.dat'])

def mv_output(path_to_output_birefCLASS, path_to_dir_save, time_crit):
    ''' ある基準時以降に作られたbirefCLASSのoutputのうち、cl_lensedのみを移動します

    Paramters:
      path_to_dir_output_birefCLASS (path) : birefCLASS の output ディレクトリ
      path_to_dir_save (path) : 移動先のディレクトリ
      datetime_crit (datetime) : 基準時刻
    '''  
    if not os.path.exists(path_to_dir_save):
        os.mkdir(path_to_dir_save)
    fos.mv_by_timecriteria(path_to_output_birefCLASS,path_to_dir_save,time_crit)
    fos.select(path_to_dir_save, ['cl_lensed.dat'])

def make_outputfilelist(path_to_dir):
    ''' 指定したディレクトリの中にあるbirefCLASSのoutputファイルについて、
        諸情報の一覧を作ります。

    Paramters:
      path_to_dir (path) : birefCLASS の output をまとめたディレクトリ
      
    returns:
      df (DataFrame) :  ファイルの一覧表
    '''
    files = fos.get_files(path_to_dir)
    filename_list = []
    beta_list = []
    mass_list = []
    neta_list = []
    for file in files:
        beta,mass,neta = get_basic_info(file)
        beta_list.append(beta)
        mass_list.append(mass)
        neta_list.append(neta)
        filename_list.append(file)
    beta_list_deg = np.rad2deg(beta_list)
    beta_list_deg_round = [round(beta, 2) for beta in beta_list_deg]
    makimaki_list = [round(beta/180) for beta in beta_list_deg]
    delta_list = [round(beta - int(beta/180)*180 -0.34, 2) for beta in beta_list_deg]
    df = pd.DataFrame({'filename':filename_list,'mass':mass_list, 'beta(rad)':beta_list, 'beta(deg)':beta_list_deg_round, 'maki':makimaki_list, 'delta(deg)':delta_list, 'neta':neta_list})
    df = df.sort_values(by=['mass', 'beta(rad)', 'neta'])   
    return df


def load_output(path_to_file):
    ''' transform CLASS output of Cl into DataFrame

    Paramters:
      path_to_file : path to the output file of CLASS

    Return:
      df_output : 
    '''
    
    df_output = pd.read_table(path_to_file,comment='#', header=None,delim_whitespace=True) 

    df_output = df_output.rename(columns={0: 'l', 1:'TT', 2:'EE',  3:'TE',4:'BB', 5:'EB', 6:'TB', 7:'phiphi', 8:'Tphi', 9:'Ephi'})
    df_output = df_output[df_output['l'] < 2501]

    lstart = df_output['l'][0]
    #TODO:l=0とl=1はTTとかも全部0でよいのか？
    for i in range(lstart):
        new_row = pd.Series([lstart-i-1, 0, 0, 0, 0, 0, 0, 0, 0, 0], index=['l','TT','EE','TE','BB','EB','TB','phiphi','Tphi','Ephi'])
        df_output = pd.concat([new_row.to_frame().transpose(), df_output], ignore_index=True)
    
    return df_output


def load_background(path_to_file):
    df_bg = pd.read_table(path_to_file,comment='#',header=None,delim_whitespace=True)
    df_bg = df_bg.rename(columns={0: 'z', 1:'proper time [Gyr]', 2:'conf.time [Mpc]',  3:'H [1/Mpc]',4:'comov.dist.', 5:'ang.diam.dist.', 6:'lum.dist.', 7:'comov.snd.hrz.', 8:'rho_g', 9:'rho_b', 10:'rho_cdm', 11:'rho_lambda', 12:'rho_ur', 13:'rho_crit', 14:'rho_tot', 15:'p_tot', 16:'p_tot_prime', 17:'gr.fac. D', 18:'gr.fac. f'})
    #ρおよびpは全て(8πG/3)が掛けられたものになっている。
    df_bg['a']=df_bg['comov.dist.']/df_bg['lum.dist.']
    df_bg['a'].fillna(1.0, inplace=True)

    return df_bg

