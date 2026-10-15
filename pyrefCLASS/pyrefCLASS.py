import os
import numpy as np
import pandas as pd
from astropy import units as u
from astropy import constants as c
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp


#//// fixed parameters ////#

# directory
#dir_pwd = '/global/homes/t/toshiyan/Work/Ongoing/rotation/code/pyrefCLASS_fix/'
#dir_ini = dir_pwd+'./initial_file/'
#dir_out = dir_pwd+'./output/'
#dir_phi = dir_pwd+'./data_axion/'
#dir_bck = dir_pwd+'./background/'

# directory
dir_ini = '/Users/salvatoresamuelesirletti/pyrefCLASS/initial_file/'
dir_out = '/Users/salvatoresamuelesirletti/pyrefCLASS/output/'
dir_phi = '/Users/salvatoresamuelesirletti/pyrefCLASS/data_axion/'
dir_bck = '/Users/salvatoresamuelesirletti/pyrefCLASS/background/'

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
path_to_birefCLASS = '/Users/salvatoresamuelesirletti/birefCLASS_distribute/birefclass/'

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

# ----------------------------------------------------------------------------------------------------------
# Background of the cosmology of the CLASS run (added in THIS version).
# The table CLASS_data_background.dat read above corresponds to ONE fixed cosmology. The conformal-time grid of the
# exported dynamics must instead come from the SAME ini that CLASS will use, because transfer.c now reads the field
# at tau = tau0_run - r, with tau0_run the conformal age of the run. Call set_background(default_ini_file) before
# exporting the dynamics (once per default ini).
import hashlib, re, shutil

_bg_ready = False
_bg_info = None

def _read_background_by_name(path):
    """Read z and the conformal time [Mpc] from a CLASS background file, locating the columns by name in the header
    (the number and order of the other columns depend on the cosmology, e.g. massive neutrinos)."""
    with open(path) as fh:
        header = [l for l in fh if l.startswith('#')]
    cz = ce = None
    for l in header:
        m1 = re.search(r'(\d+):z\b', l)
        m2 = re.search(r'(\d+):conf\. ?time', l)
        if m1: cz = int(m1.group(1)) - 1
        if m2: ce = int(m2.group(1)) - 1
    assert cz is not None and ce is not None, 'columns z / conf. time not found in the header of ' + path
    data = np.loadtxt(path)
    return data[:, cz], data[:, ce]

def set_background(default_ini_file, tau_reio=0.054, ln10As=3.0445, Neta=10000, force=False):
    """Build (or load from the cache in dir_bck) the background table of the cosmology of default_ini_file by running
    CLASS once, and use it for the conversion z <-> conformal time in the export functions."""
    global z2eta, eta2z, _bg_ready, _bg_info
    _bg_ready, _bg_info = False, None
    with open(dir_ini + default_ini_file, 'rb') as fh:
        ini_hash = hashlib.md5(fh.read()).hexdigest()[:8]
    cache = dir_bck + 'bg_' + default_ini_file.replace('.ini', '') + '_' + ini_hash + '.dat'
    if force or not os.path.exists(cache):
        tmp = dir_bck + 'tmp_run/'
        os.makedirs(tmp, exist_ok=True)
        tag = 'bgtmp_' + default_ini_file.replace('.ini', '')
        dummy = dir_phi + tag
        np.savetxt(dummy, np.zeros(Neta))        # dummy dynamics file (alpha_0 = 0): CLASS opens it anyway
        create_inifile(0., tag, Neta, tag, default_ini_file=default_ini_file, tau_reio=tau_reio, ln10As=ln10As,
                       outroot=tmp, phi_dynamics_dir=dir_phi)
        os.system(path_to_birefCLASS + '/class ' + dir_ini + tag + '.ini ' + path_to_birefCLASS + '/cl_ref_biref.pre > /dev/null')
        bg_out = tmp + tag + '_background.dat'
        assert os.path.exists(bg_out), 'CLASS did not write ' + bg_out
        shutil.copyfile(bg_out, cache)
        os.system('rm -f ' + tmp + tag + '_* ' + dir_ini + tag + '.ini ' + dummy)
    z, eta = _read_background_by_name(cache)
    z2eta = interp1d(z, eta, kind='linear', fill_value='extrapolate')
    eta2z = interp1d(eta, z, kind='linear', fill_value='extrapolate')
    _bg_ready, _bg_info = True, (default_ini_file, cache)
    print('pyrefCLASS background set from', default_ini_file, '-> eta_0 =', float(z2eta(0)), 'Mpc')

def _check_background():
    if not _bg_ready:
        raise RuntimeError("Call set_background(default_ini_file) before exporting the dynamics: the conformal-time grid "
                           "must come from the same ini as the CLASS run.")

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

    _check_background()
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

    _check_background()
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
    


# ----------------------------------------------------------------------------------------------------------
# Safety check (added in THIS version): the default ini used by create_inifile must be the one of set_background.
import functools
_create_inifile_orig = create_inifile

@functools.wraps(_create_inifile_orig)
def create_inifile(*args, **kwargs):
    ini = kwargs.get('default_ini_file', args[4] if len(args) > 4 else None)
    if ini is not None and _bg_info is not None and ini != _bg_info[0]:
        raise RuntimeError("create_inifile uses " + ini + " but the background was set for " + _bg_info[0] +
                           ": call set_background('" + ini + "') first.")
    return _create_inifile_orig(*args, **kwargs)
