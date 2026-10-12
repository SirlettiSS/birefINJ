import numpy as np
import os
import sys
import uuid
import shutil
import itertools
from astropy.io import fits
from multiprocessing import Pool, cpu_count

# ==============================================================================
# GRID SEARCH SCRIPT FOR BIREFCLASS: HPC ARCHITECTURE (EE-ONLY LIKELIHOOD)
# ==============================================================================
# FROM MCMC TO GRID SEARCH:
# The MCMC approach (emcee) failed to properly sample the highly degenerate 
# parameter space in the EE-only channel, leading to statistical artifacts and 
# walkers clustering against prior boundaries. This script replaces the stochastic 
# sampling with a deterministic 2D Grid Search (Profile Likelihood) to map the 
# exact chi-square topology without sampling noise.
#
# [HPC-Ready Approach]
# 1. MULTIPROCESSING: Grid points are evaluated in parallel across available CPU cores.
# 2. RAM DISK (/dev/shm): All temporary files (.ini, .dat, CLASS outputs) are written 
#    to the node's volatile RAM and not on the physical disk. This prevents severe 
#    I/O bottlenecks and file system lock-ups on the cluster.
# 3. WORKER ISOLATION: Each grid point evaluation generates a unique UUID folder 
#    in RAM to prevent any race conditions or file overwriting between concurrent cores.
# 4. TEMPLATE REDIRECTION: We temporarily redirect `rot.dir_ini` to our RAM folder 
#    so pyrefCLASS writes the .ini there, without needing to alter the core library.
# 5. SMART CORES: Automatically detects PBS allocation to avoid CPU overloading.
#
# [2D GRID VERSION]
# Mass (log10_ma) and injection window (z_start, z_end) are fixed via environment variables.
# The grid systematically explores the (xi, tau) plane.
#
# [LIKELIHOOD]
# The code evaluates the exact chi-square using a Gaussian likelihood in Dl_EE, 
# built from the Planck NPIPE lollipop data (low-ell EE spectrum, lmax=30) and 
# its associated covariance matrix. This EE-only run mainly serves as a diagnostic 
# baseline to expose the mathematical degeneracy before extending to EE+EB+BB.
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
outroot_grid    = '/home/sirlettiss/biref_from_injection/biref_iso_inject/EE_grid_search/grid_search_output/'
RAM_DISK_DIR    = '/dev/shm'

os.makedirs(outroot_grid, exist_ok=True)

sys.path.append(dir_pyrefclass)
sys.path.append('/home/sirlettiss/biref_from_injection/biref_iso_inject/')
import pyrefCLASS as rot
import local

ORIGINAL_INI_DIR = rot.dir_ini

# ==============================================================================
# 2. OBSERVATIONAL DATA LOADING (EE ONLY)
# ==============================================================================
def data_npipe(Dl=False, lmax=30, spec='all'):
    lb   = np.linspace(2, lmax, lmax - 1)
    nbin = lmax - 1
    lfac = lb * (lb + 1) / 2 / np.pi if Dl else np.ones(nbin)

    cldat = np.loadtxt(os.path.join(data_path, 'cl_lolEB_NPIPE.dat'), unpack=True)[1:4, :nbin] * lfac[None, :]
    clcov = fits.getdata(os.path.join(data_path, 'clcov_lolEB_NPIPE.fits'))

    cl = cldat[0]
    cov = clcov[0:nbin, 0:nbin]

    if Dl:
        cov = cov * np.outer(lfac, lfac)

    vcl  = cov.diagonal()
    icov = np.linalg.inv(cov)
    scl  = np.sqrt(vcl)

    return lb, cl, scl, icov, vcl

lb, Dl_EE_obs, scl_EE, icov_EE, vcl_EE = data_npipe(Dl=True, lmax=30, spec='EE')

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

        # Exact Chi-Square Computation
        lcl = np.loadtxt(ofile, unpack=True, max_rows=3000)
        nbin = len(Dl_EE_obs)
        
        # Apply thermodynamic scaling
        EE_Dl_th = lcl[2][:nbin] * TCMB_muK**2
        
        delta = EE_Dl_th - Dl_EE_obs
        chi2 = delta @ icov_EE @ delta
        
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

    print(f"=== STARTING 2D GRID SEARCH (EE-ONLY) ===")
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

    run_tag = f'Grid2D_EE_logma{LOG10_MA}_z{int(Z_START)}-{int(Z_END)}_{N_POINTS_XI}x{N_POINTS_TAU}'
    out_file = os.path.join(outroot_grid, f'{run_tag}.npz')
    
    np.savez(out_file, xi=xi_grid, tau=tau_grid, chi2=chi2_grid)
    
    print(f"=== GRID SEARCH COMPLETED ===")
    print(f"Results successfully saved in: {out_file}")