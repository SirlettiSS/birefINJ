import numpy as np
import os
import sys
import uuid
import shutil
import itertools
from astropy.io import fits
from multiprocessing import Pool, cpu_count

# ==============================================================================
# GRID SEARCH SCRIPT FOR BIREFCLASS: HPC ARCHITECTURE (JOINT LIKELIHOOD)
# ==============================================================================
# [JOINT LIKELIHOOD: EE + EB + BB]
# This script evaluates the exact chi-square topology over a 2D (xi, tau) grid.
# Unlike the EE-only run which suffers from severe mathematical degeneracy, 
# this version uses the joint polarization likelihood (EE, EB, BB) from the 
# Planck NPIPE lollipop data. 
# 
# The energy removed from the EE channel due to the axion coupling is transferred 
# into the EB and BB channels. By fitting all three simultaneously, the parameter 
# space becomes properly constrained.
#
# [HPC-Ready Approach]
# 1. MULTIPROCESSING: Grid points are evaluated in parallel across CPU cores.
# 2. RAM DISK (/dev/shm): All temporary files (.ini, .dat, CLASS outputs) are written 
#    to volatile RAM, preventing I/O bottlenecks.
# 3. WORKER ISOLATION: Unique UUID folders prevent race conditions.
# ==============================================================================

# ==============================================================================
# 0. ENVIRONMENT & GRID PARAMETERS (Read dynamically from PBS script)
# ==============================================================================
try:
    Z_START  = float(os.environ.get('Z_START'))
    Z_END    = float(os.environ.get('Z_END'))
    XI_MAX   = float(os.environ.get('XI_MAX'))
    LOG10_MA = float(os.environ.get('LOG10_MA', -30.0))
except TypeError:
    print("CRITICAL ERROR: Z_START, Z_END, or XI_MAX not defined in the environment.")
    sys.exit(1)

# Grid settings (with fallbacks if not provided)
N_POINTS_XI  = int(os.environ.get('N_POINTS_XI', 40))
N_POINTS_TAU = int(os.environ.get('N_POINTS_TAU', 40))
TAU_MIN      = float(os.environ.get('TAU_MIN', 0.04))
TAU_MAX      = float(os.environ.get('TAU_MAX', 0.12))

# ==============================================================================
# 1. PATH SETUP (IPMU Cluster)
# ==============================================================================
dir_pyrefclass  = '/home/sirlettiss/biref_from_injection/pyrefCLASS'
path_birefclass = '/home/sirlettiss/biref_from_injection/birefclass'
data_path       = '/home/sirlettiss/biref_from_injection/lollipop_data/planck_2020/lollipop/'
outroot_grid    = '/home/sirlettiss/biref_from_injection/biref_iso_inject/joint_grid_search/grid_search_output/'
RAM_DISK_DIR    = '/dev/shm'

os.makedirs(outroot_grid, exist_ok=True)

sys.path.append(dir_pyrefclass)
sys.path.append('/home/sirlettiss/biref_from_injection/biref_iso_inject/')
import pyrefCLASS as rot
import local

ORIGINAL_INI_DIR = rot.dir_ini

# ==============================================================================
# 2. OBSERVATIONAL DATA LOADING (JOINT EE + EB + BB)
# ==============================================================================
def data_npipe(Dl=False, lmax=30, spec='all'):
    lb   = np.linspace(2, lmax, lmax - 1)
    nbin = lmax - 1
    lfac = lb * (lb + 1) / 2 / np.pi if Dl else np.ones(nbin)

    cldat = np.loadtxt(os.path.join(data_path, 'cl_lolEB_NPIPE.dat'), unpack=True)[1:4, :nbin] * lfac[None, :]
    clcov = fits.getdata(os.path.join(data_path, 'clcov_lolEB_NPIPE.fits'))
    nell  = len(clcov) // 3

    if spec == 'EE':
        cl = cldat[0];  cov = clcov[0:nbin, 0:nbin]
    elif spec == 'BB':
        cl = cldat[1];  cov = clcov[nell:nell+nbin, nell:nell+nbin]
    elif spec == 'EB':
        cl = cldat[2];  cov = clcov[2*nell:2*nell+nbin, 2*nell:2*nell+nbin]
    elif spec == 'all':
        cl   = cldat.copy()
        cbcov = np.zeros((3*nbin, 3*nbin))
        for t1 in range(3):
            for t2 in range(3):
                cbcov[t1*nbin:(t1+1)*nbin, t2*nbin:(t2+1)*nbin] = \
                    clcov[t1*nell:t1*nell+nbin, t2*nell:t2*nell+nbin]
        cov = cbcov
    else:
        raise ValueError(f"spec='{spec}' not valid.")

    if Dl:
        if spec == 'all':
            lfac_all = np.tile(lfac, 3)
            cov = cov * np.outer(lfac_all, lfac_all)
        else:
            cov = cov * np.outer(lfac, lfac)

    vcl  = cov.diagonal()
    icov = np.linalg.inv(cov)
    
    if spec == 'all':
        scl = np.sqrt(vcl.reshape(-1, nbin))
    else:
        scl = np.sqrt(vcl)

    return lb, cl, scl, icov, vcl

# Load all three spectra and flatten the observation matrix
lb, Dl_all_obs, scl_all, icov_joint, vcl_joint = data_npipe(Dl=True, lmax=30, spec='all')
Dl_obs_vector = Dl_all_obs.flatten()

# ==============================================================================
# 3. SINGLE POINT EVALUATION (Isolated Topology in RAM disk)
# ==============================================================================
Neta         = 10000
beta0        = 0.34 * np.pi / 180.
n, tf, eps   = 1, 100, 1e-8
theta_ini    = 0.01
bphi0_ini    = theta_ini * tf
TCMB_muK     = 2.7255e6
const_planck = 2.1e-9 * np.exp(-2 * 0.054)

def compute_grid_point(params):
    xi, tau = params
    ln10As = np.log(1e10 * const_planck * np.exp(2 * tau))
    
    worker_id = uuid.uuid4().hex[:8]
    base_tag  = f'grid_{worker_id}'
    worker_ram_dir = os.path.join(RAM_DISK_DIR, f'grid_{worker_id}') + os.sep
    
    os.makedirs(worker_ram_dir, exist_ok=True)
    template_name = 'biref_default_inject.ini'
    shutil.copy(os.path.join(ORIGINAL_INI_DIR, template_name), worker_ram_dir)
    
    rot.dir_ini = worker_ram_dir
    
    try:
        dyn_file = f"{base_tag}_phi.dat"
        dyn_path = os.path.join(worker_ram_dir, dyn_file)
        ofile    = os.path.join(worker_ram_dir, base_tag + '_cl.dat')

        rot.export_axion_dynamics_a(Neta, local.EoM_phi,
                                    [LOG10_MA, tf, n, xi, eps, Z_START, Z_END, bphi0_ini], dyn_path)

        rot.create_inifile(beta0, base_tag, Neta, dyn_file,
                           default_ini_file=template_name,
                           tau_reio=tau, ln10As=ln10As,
                           outroot=os.path.join(worker_ram_dir, ''),
                           phi_dynamics_dir=worker_ram_dir)

        ini_file_path = os.path.join(worker_ram_dir, base_tag + '.ini')
        os.system(f"{path_birefclass}/class {ini_file_path} "
                  f"{path_birefclass}/cl_ref_biref.pre > /dev/null 2>&1")

        # Exact Joint Chi-Square Computation
        lcl = np.loadtxt(ofile, unpack=True, max_rows=3000)
        nbin = len(Dl_obs_vector) // 3
        
        # Apply thermodynamic scaling to all three channels
        Dl_EE_th = lcl[2][:nbin] * TCMB_muK**2
        Dl_BB_th = lcl[4][:nbin] * TCMB_muK**2
        Dl_EB_th = lcl[5][:nbin] * TCMB_muK**2
        
        # Concatenate the theoretical vectors in the exact order of the covariance matrix
        Dl_th_vector = np.concatenate([Dl_EE_th, Dl_BB_th, Dl_EB_th])
        
        delta = Dl_th_vector - Dl_obs_vector
        chi2 = delta @ icov_joint @ delta
        
        return xi, tau, chi2

    except Exception as e:
        print(f"Warning: CLASS failed for xi={xi:.5e}, tau={tau:.4f}. Error: {e}")
        return xi, tau, np.inf
        
    finally:
        rot.dir_ini = ORIGINAL_INI_DIR
        shutil.rmtree(worker_ram_dir, ignore_errors=True)

# ==============================================================================
# 4. GRID SEARCH ENGINE
# ==============================================================================
if __name__ == "__main__":
    
    # 1. Grid Definition
    xi_vals  = np.linspace(0.0, XI_MAX, N_POINTS_XI)
    tau_vals = np.linspace(TAU_MIN, TAU_MAX, N_POINTS_TAU)
    
    # Create all possible (xi, tau) combinations
    grid_points = list(itertools.product(xi_vals, tau_vals))

    # 2. Core Management
    pbs_cpus = os.environ.get('PBS_CPUS')
    available_cpus = int(pbs_cpus) if pbs_cpus else cpu_count()

    print(f"=== STARTING 2D GRID SEARCH (JOINT EE+EB+BB) ===")
    print(f"Log10_ma: {LOG10_MA}")
    print(f"Redshift Window: z={Z_START} -> z={Z_END}")
    print(f"Xi Grid:  [0.0, {XI_MAX}] ({N_POINTS_XI} points)")
    print(f"Tau Grid: [{TAU_MIN}, {TAU_MAX}] ({N_POINTS_TAU} points)")
    print(f"Total evaluations: {len(grid_points)}")
    print(f"Allocated cores: {available_cpus}")

    # 3. Parallel Execution
    with Pool(processes=available_cpus) as pool:
        results = pool.map(compute_grid_point, grid_points)

    # 4. Data Restructuring and Saving
    results_array = np.array(results)
    
    # Reshape linear results into 2D matrices for easy contour plotting
    xi_grid   = results_array[:, 0].reshape(N_POINTS_XI, N_POINTS_TAU)
    tau_grid  = results_array[:, 1].reshape(N_POINTS_XI, N_POINTS_TAU)
    chi2_grid = results_array[:, 2].reshape(N_POINTS_XI, N_POINTS_TAU)

    # Dynamic run_tag to specify JOINT analysis
    run_tag = f'Grid2D_JOINT_logma{LOG10_MA}_z{int(Z_START)}-{int(Z_END)}_{N_POINTS_XI}x{N_POINTS_TAU}'
    out_file = os.path.join(outroot_grid, f'{run_tag}.npz')
    
    np.savez(out_file, xi=xi_grid, tau=tau_grid, chi2=chi2_grid)
    
    print(f"=== GRID SEARCH COMPLETED ===")
    print(f"Results successfully saved in: {out_file}")