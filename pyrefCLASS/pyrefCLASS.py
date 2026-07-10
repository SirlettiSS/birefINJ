import os
import numpy as np
import pandas as pd
from astropy import units as u
from astropy import constants as c
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp


#//// fixed parameters ////#

# directory

dir_ini = '/leonardo_work/INF26_indark/ssirlett/biref_from_injection/pyrefCLASS/initial_file/'
dir_out = '/leonardo_work/INF26_indark/ssirlett/biref_from_injection/pyrefCLASS/output/'
dir_phi = '/leonardo_work/INF26_indark/ssirlett/biref_from_injection/pyrefCLASS/data_axion/'
dir_bck = '/leonardo_work/INF26_indark/ssirlett/biref_from_injection/pyrefCLASS/background/'

if not os.path.exists(dir_ini):
    os.makedirs(dir_ini)
else:
    pass
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
else:
    pass
if not os.path.exists(dir_phi):
    os.makedirs(dir_phi)
else:
    pass

# birefclass
path_to_birefCLASS = '/leonardo_work/INF26_indark/ssirlett/biref_from_injection/birefclass'

# constant
Gyr2Mpc = u.Glyr.to(u.Mpc)
eV2MpcInv = (1*u.eV/(c.hbar*c.c)).to(1/u.Mpc).value #確認すること

df_bg = pd.read_table(dir_bck+'CLASS_data_background.dat',comment='#',header=None,sep=r'\s+')
df_bg = df_bg.rename(columns={0: 'z', 1:'proper time [Gyr]', 2:'conf.time [Mpc]',  3:'H [1/Mpc]',4:'comov.dist.', 5:'ang.diam.dist.', 6:'lum.dist.', 7:'comov.snd.hrz.', 8:'rho_g', 9:'rho_b', 10:'rho_cdm', 11:'rho_lambda', 12:'rho_ur', 13:'rho_crit', 14:'rho_tot', 15:'p_tot', 16:'p_tot_prime', 17:'gr.fac. D', 18:'gr.fac. f'})
#ρおよびpは全て(8πG/3)が掛けられたものになっている。
df_bg['a']=df_bg['comov.dist.']/df_bg['lum.dist.']
#df_bg['a'].fillna(1.0, inplace=True)
df_bg['a'] = df_bg['a'].fillna(1.0)

z2eta = interp1d(df_bg['z'], df_bg['conf.time [Mpc]'], kind='linear', fill_value='extrapolate')
eta2z = interp1d(df_bg['conf.time [Mpc]'], df_bg['z'], kind='linear', fill_value='extrapolate')
z2time = interp1d(df_bg['z'], Gyr2Mpc*df_bg['proper time [Gyr]'])
#time2z = interp1d(Gyr2Mpc*df_bg['proper time [Gyr]'], df_bg['z'])
eta2time = interp1d(df_bg['conf.time [Mpc]'], Gyr2Mpc*df_bg['proper time [Gyr]']) #確認すること
time2Hubble = interp1d( Gyr2Mpc*df_bg['proper time [Gyr]'], df_bg['H [1/Mpc]'], kind='linear', fill_value='extrapolate')

    

def EoM_phi_mass_potential(t, y, mass):
    '''
    EoM equation of mass potential model
    \ddot{\phi} + 3H(t)\dot{\phi}+ m^2\phi = 0

    Args:
    t (float) : proper time 
    y [float, float] : [0th derivative , 1st derivative]
    mass (float) : the mass of ALP (eV)

    Returns:
    [float, float] :  [1st derivative, 2nd derivative] 
    '''

    return [ y[1], - 3 * time2Hubble(t) * y[1] - (eV2MpcInv**2) * (mass**2) * y[0]] 


def EoM_phi_cos_potential(t, y, mass, f):
    '''
    EoM equation of Axion (cos) potential model
    \ddot{\phi} + 3H(t)\dot{\phi}+ m^2*f*\sin(\phi/f)= 0

    Args:
    t (float) : proper time 
    y [float, float] : [0th derivative , 1st derivative]
    mass (float) : the mass of ALP (eV)
    f (float) : the constant of the potential

    Returns:
    [float, float] :  [1st derivative, 2nd derivative] 
    '''
        
    return [ y[1], - 3 * time2Hubble(t) * y[1] - (eV2MpcInv**2) * (mass**2) * f * np.sin(f * y[0])] 


def filename_axion_dynamics(Neta,args):
    
    # define filename
    paramname = ''
    for i, arg in enumerate(args):
        name = '_p' + str(i) + '_' + str(arg)
        paramname += name
    
    return 'axionDynamics' + paramname + '_N' + str(Neta) + '.dat'
    

def export_axion_dynamics_a(Neta, model, args, output_filename, initial_redshift=10**6, phi0ini=1, phi1ini=0, atol = 1e-10, rtol = 1e-5):
    '''
    アクシオンのdynamicsを./data_axion/以下にdatファイルで出力します

    Args:
    etalist (list) : conformal time list
    zlist (list) : redshift list
    model (func) : EoM as a function of scale factor
    args (list) : 上記 EoM に必要なパラメータのリスト
    initial_redshift (int or float) : 赤方偏移で表した初期時刻
    phi0ini (int or float) : アクシオン場の初期値
    phi1ini (int or float) : アクシオン場の1階微分の初期値
    atol (float)
    rtol (float)
    ow (bool): overwrite or not

    Returns:
    filename (str): ***.dat 
    '''

    eta0 = z2eta(0) # conformal time at present
    etalist = eta0 * np.arange(1, Neta + 1) / Neta
    alist = 1./(1.+eta2z(etalist))

    if not os.path.exists(output_filename):
        
        # time steps
        aini = 1./(1.+initial_redshift)
        atoday = 1.

        # solve EoM
        sol = solve_ivp(model, [aini,atoday], [phi0ini,phi1ini], t_eval=alist, args=args, atol=atol, rtol=rtol)
        solphi0 = sol['y'][0]
        phi0list = (solphi0 - solphi0[-1])/(phi0ini - solphi0[-1]) #要確認
    
        # save to file
        np.savetxt(output_filename, phi0list)


def export_axion_dynamics(Neta, model, args, initial_redshift=10**6, phi0ini=1, phi1ini=0, atol = 1e-10, rtol = 1e-5):
    '''
    アクシオンのdynamicsを./data_axion/以下にdatファイルで出力します

    Args:
    timelist (list) : the time list of proper time 
    model (func) : EoM の微分方程式の関数
    args (list) : 上記 EoM に必要なパラメータのリスト
    initial_redshift (int or float) : 赤方偏移で表した初期時刻
    phi0ini (int or float) : アクシオン場の初期値
    phi1ini (int or float) : アクシオン場の1階微分の初期値
    atol (float)
    rtol (float)
    ow (bool): overwrite or not

    Returns:
    filename (str): ***.dat 
    '''

    eta0 = z2eta(0) # conformal time at present
    etalist = eta0 * np.arange(1, Neta + 1) / Neta
    timelist = eta2time(etalist)

    filename = filename_axion_dynamics(Neta,args)

    if not os.path.exists(dir_phi+filename):
        
        # time steps
        timeini = z2time(initial_redshift)
        timetoday = z2time(0)

        # solve EoM
        sol = solve_ivp(model, [timeini,timetoday], [phi0ini,phi1ini], t_eval=timelist, args=args, atol=atol, rtol=rtol)
        solphi0 = sol['y'][0]
        phi0list = (solphi0 - solphi0[-1])/(phi0ini - solphi0[-1]) #要確認
    
        # save to file
        np.savetxt(dir_phi+filename, phi0list)


def create_inifile(alpha0, base_tag, Neta, filename_dynamics, default_ini_file='biref_default_makimaki.ini', tau_reio = None, ln10As=None,
                   outroot = None, phi_dynamics_dir = None,
                  ):
    '''
    BirefCLASSのinitial fileを作ります

    Args:
    alpha0 (float) : birefringence angle [rad]
    Neta (int) : 時刻の刻み数
    filename (str) : export_axion_dynamics で出力されたファイル名 (***.dat)

    Returns:
    (str): 出力されるinitial file のファイル名(***.ini)
    '''

    inifilename = base_tag + '.ini'
    cwd = os.getcwd()

    if outroot is None:
        oroot = cwd + '/output/' + base_tag
    else:
        oroot = outroot + base_tag

    # 1. 元の内容を読み込む & tau_reio の行を除去
    with open(dir_ini + default_ini_file, 'r') as f:
        lines = f.readlines()
    # tau_reio の行を除外
    if tau_reio is not None:
        lines = [line for line in lines if not line.strip().startswith("tau_reio")]
    
    if ln10As is not None:
        lines = [l for l in lines if not l.strip().startswith("ln10^{10}A_s")]

    # axion dynamics のファイルパスを決定
    if phi_dynamics_dir is None:
        dynamics_path = os.path.join(cwd, dir_phi, filename_dynamics)
    else:
        dynamics_path = os.path.join(phi_dynamics_dir, filename_dynamics)

    # iniファイルのフルパスを作成
    ini_path = os.path.join(dir_ini, inifilename)

    # iniファイルに必要な設定を出力
    with open(ini_path, 'w') as inifile:
        config_text = (
            f"filename_axion_dynamics={dynamics_path}\n"
            f"length_resolution_axion_dynamics={Neta}\n"
            f"isotropic_alpha_0={alpha0}\n"
            f"is_m_minus=no\n"
            f"root={oroot+'_'}\n"
        )
        inifile.write(config_text)
    

        # 残しておく元の行を追記
        for line in lines:
            inifile.write(line)

        # tau_reio の値があれば最後に追加
        if tau_reio is not None:
            inifile.write(f"tau_reio={float(tau_reio)}\n")
            
        if ln10As is not None:                              
            inifile.write(f"ln10^{{10}}A_s = {float(ln10As)}\n")

    '''
    with open(default_ini, 'r') as default_file:
        default_inifile = default_file.read()

    with open(dir_ini+inifilename, 'w') as inifile:
        inifile.write(
            f"filename_axion_dynamics={cwd}/{dir_phi}/{filename_dynamics}\n" #絶対パスで指定するようにした
            f"length_resolution_axion_dynamics={Neta}\n"
            f"isotropic_alpha_0={alpha0}\n"
            f"is_m_minus=no\n"
            f"root={cwd}/output/{base_tag}\n" #絶対パスで指定するようにした
            f"{default_inifile}"            
        )
        if tau_reio is not None:
            inifile.write(f"tau_reio={float(tau_reio)}\n")
    '''
    


