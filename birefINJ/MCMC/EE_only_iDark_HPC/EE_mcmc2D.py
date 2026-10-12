import numpy as np
import os
import sys
import uuid
import shutil
from astropy.io import fits
import emcee
import h5py
from multiprocessing import Pool, cpu_count

# ==============================================================================
# MCMC SCRIPT FOR BIREFCLASS: LOCAL VS HPC ARCHITECTURE
# ==============================================================================

# [2D VERSION]
# Mass log10_ma is fixed
LOG10_MA_FIXED = -30.5

# Read environment variables passed by SLURM (HPC Architecture)
try:
    z_start = float(os.environ.get('Z_START'))
    z_end   = float(os.environ.get('Z_END'))
    xi_max  = float(os.environ.get('XI_MAX'))
except TypeError:
    print("ERROR: Environment variables Z_START, Z_END, or XI_MAX not found.")
    sys.exit(1)

# ==============================================================================
# 0. ENVIRONMENT SETTINGS (Local vs HPC)
# ==============================================================================

RUN_ON_HPC = True

if RUN_ON_HPC:
    dir_pyrefclass  = '/home/sirlettiss/biref_from_injection/pyrefCLASS'
    path_birefclass = '/home/sirlettiss/biref_from_injection/birefclass'
    data_path       = '/home/sirlettiss/biref_from_injection/lollipop_data/planck_2020/lollipop/'
    outroot_mcmc    = '/home/sirlettiss/biref_from_injection/biref_iso_inject/logma_30.5_diffredbins/output/'
    os.makedirs(outroot_mcmc, exist_ok=True)
    RAM_DISK_DIR    = '/dev/shm'
else:
    dir_pyrefclass  = '/Users/salvatoresamuelesirletti/pyrefCLASS'
    path_birefclass = '/Users/salvatoresamuelesirletti/birefCLASS_distribute/birefclass'
    data_path       = '/Users/salvatoresamuelesirletti/lollipop_data/planck_2020/lollipop/'
    outroot_mcmc    = '/Users/salvatoresamuelesirletti/Desktop/Work/Biref_iso_inject/output/'
    import tempfile
    RAM_DISK_DIR = tempfile.gettempdir()

sys.path.append(dir_pyrefclass)
sys.path.append('/home/sirlettiss/biref_from_injection/biref_iso_inject/')
import pyrefCLASS as rot
import local

ORIGINAL_INI_DIR = rot.dir_ini

# ==============================================================================
# 1. DATA LOADING (EE ONLY)
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

    # From Cov in Cl to Cov in Dl (lfac * lfac^T)
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

# Load only the EE spectrum
lb, Dl_EE_obs, scl_EE, icov_EE, vcl_EE = data_npipe(Dl=True, lmax=30, spec='EE')

# ==============================================================================
# 2. CLASS COMPUTATION (thread-safe)
# ==============================================================================

Neta         = 10000
beta0        = 0.34 * np.pi / 180.
n, tf, eps   = 1, 100, 1e-8
theta_ini    = 0.01
bphi0_ini    = theta_ini * tf
TCMB_muK     = 2.7255e6
const_planck = 2.1e-9 * np.exp(-2 * 0.054)

def compute_cls(xi, log10_ma, tau):
    ln10As    = np.log(1e10 * const_planck * np.exp(2 * tau))
    worker_id = uuid.uuid4().hex[:8]
    base_tag  = f'inj_{worker_id}'

    worker_ram_dir = os.path.join(RAM_DISK_DIR, f'mcmc_{worker_id}') + os.sep
    os.makedirs(worker_ram_dir, exist_ok=True)

    template_name = 'biref_default_inject.ini'
    shutil.copy(os.path.join(ORIGINAL_INI_DIR, template_name), worker_ram_dir)

    rot.dir_ini = worker_ram_dir

    try:
        dyn_file = f"{base_tag}_phi.dat"
        dyn_path = os.path.join(worker_ram_dir, dyn_file)
        ofile    = os.path.join(worker_ram_dir, base_tag + '_cl.dat')

        # 8-element array including z_start, z_end (from SLURM) and bphi0_ini
        rot.export_axion_dynamics_a(Neta, local.EoM_phi,
                                    [log10_ma, tf, n, xi, eps, z_start, z_end, bphi0_ini], dyn_path)

        rot.create_inifile(beta0, base_tag, Neta, dyn_file,
                           default_ini_file=template_name,
                           tau_reio=tau, ln10As=ln10As,
                           outroot=os.path.join(worker_ram_dir, ''),
                           phi_dynamics_dir=worker_ram_dir)

        ini_file_path = os.path.join(worker_ram_dir, base_tag + '.ini')

        os.system(f"{path_birefclass}/class {ini_file_path} "
                  f"{path_birefclass}/cl_ref_biref.pre > /dev/null 2>&1")

        lcl = np.loadtxt(ofile, unpack=True, max_rows=3000)
        lcl[1:] *= TCMB_muK**2

    finally:
        rot.dir_ini = ORIGINAL_INI_DIR
        shutil.rmtree(worker_ram_dir, ignore_errors=True)

    return lcl

# ==============================================================================
# 3. LIKELIHOOD (2D: xi and tau only, EE spectrum)
# ==============================================================================

def log_prior(params):
    xi, tau = params
    # Modified using xi_max read from SLURM
    if (0 < xi < xi_max) and (0.04 < tau < 0.12):
        return 0.0
    return -np.inf

def log_probability(params, Dl_obs, icov):
    lp = log_prior(params)
    if not np.isfinite(lp):
        return -np.inf
    try:
        lcl   = compute_cls(params[0], LOG10_MA_FIXED, params[1])
        nbin  = len(Dl_obs)
        # Using lcl[2] which corresponds to EE
        delta = lcl[2][:nbin] - Dl_obs
        return lp - 0.5 * (delta @ icov @ delta)
    except Exception as e:
        print(f"Warning: CLASS failed for params {params}: {e}")
        return -np.inf

# ==============================================================================
# 4. MCMC with HDF5 backend (automatic resume)
# ==============================================================================

if __name__ == "__main__":

    ndim      = 2
    nwalkers  = 32
    max_steps = 5000

    os.makedirs(outroot_mcmc, exist_ok=True)

    # Dynamic file saving nomenclature for EE
    run_tag          = f'2D_EE_logma{LOG10_MA_FIXED}_z{int(z_start)}-{int(z_end)}_ximax{xi_max:.2e}'
    backend_filename = os.path.join(outroot_mcmc, f'biref_mcmc_chain_{run_tag}.h5')
    backend          = emcee.backends.HDFBackend(backend_filename)

    if os.path.exists(backend_filename) and backend.iteration > 0:
        print(f"Resuming from step {backend.iteration}...")
        initial_pos = None
        steps_to_do = max_steps - backend.iteration
    else:
        print("Starting new chain...")
        backend.reset(nwalkers, ndim)
        
        # Lowered initial guess optimized for smaller xi_max limits
        p0          = np.array([1e-3, 0.08]) 
        initial_pos = p0 + 1e-4 * np.random.randn(nwalkers, ndim)
        steps_to_do = max_steps

    slurm_cpus     = os.environ.get('SLURM_CPUS_PER_TASK')
    available_cpus = int(slurm_cpus) if slurm_cpus else cpu_count()
    n_cores        = min(nwalkers, available_cpus)

    print(f"[2D MCMC EE ONLY] log10_ma = {LOG10_MA_FIXED} fixed")
    print(f"Starting MCMC: {n_cores} cores, {nwalkers} walkers, {steps_to_do} steps to go.")

    with Pool(processes=n_cores) as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability,
                                        args=(Dl_EE_obs, icov_EE),
                                        pool=pool,
                                        backend=backend)
        if steps_to_do > 0:
            sampler.run_mcmc(initial_pos, steps_to_do, progress=True)
        else:
            print("Chain already complete!")

    print("MCMC finished!")
    print(f"Chain shape:              {sampler.get_chain().shape}")
    print(f"Mean acceptance fraction: {np.mean(sampler.acceptance_fraction):.3f}")

    np.save(os.path.join(outroot_mcmc, f'chain_{run_tag}_{max_steps}steps.npy'),    sampler.get_chain())
    np.save(os.path.join(outroot_mcmc, f'log_prob_{run_tag}_{max_steps}steps.npy'), sampler.get_log_prob())

    print(f"Chain saved: chain_{run_tag}_{max_steps}steps.npy")