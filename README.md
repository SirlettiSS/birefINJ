# Cosmic Birefringence from Axion Energy Injection

This module implements a model for isotropic cosmic birefringence in which
an energy transfer from matter to an axion-like particle (ALP) occurs during
a specific period after recombination. The resulting evolution of the axion
field induces a rotation of the CMB polarization plane, transferring power
from EE/BB into the EB spectrum.


## Requirements

To run this module, you need to install:

- [pyrefCLASS](https://github.com/deer-astron/pyrefCLASS)
- [birefCLASS](https://github.com/toshiyan/ongoing/tree/master/rotation/code/birefclass_lens)

Both folders are already included in this repository. If you run into issues
compiling birefCLASS on your HPC cluster, `recompile_CLASS.sh` should help.

The data used for the MCMC are the Planck 2020 NPIPE lollipop likelihood data:
[https://github.com/planck-npipe/lollipop](https://github.com/planck-npipe/lollipop)

Before running any script, create a directory to store the output files and
update the module/data paths accordingly in each script.


## Repository structure

### `birefINJ/`

- `local.py`: defines the injection dynamics (axion equation of motion).
- Additional Jupyter notebooks for various analyses (kick-out limit, C_ell
  with injection vs. NPIPE data, kick-out limit across different redshift
  bins, etc.).
- `output/`: contains the C_ell spectra computed through birefCLASS.

### `MCMC/`

Contains the scripts to run the 2D MCMC, either using the EE spectrum only
or jointly with EE+BB+EB. Results are stored in `MCMC/outputs/`, where a
dedicated notebook (`result_analysis/`) helps read and visualize the chains.


### `grid_search/`

Contains the scripts to run a 2D grid-based profile likelihood over (ξ, τ) at
fixed axion mass, either using the EE spectrum only or jointly with EE+BB+EB.
The exact χ² is evaluated on each grid point (parallelized, with temporary
CLASS files written to `/dev/shm`), which avoids the sampling issues of the
MCMC in the degenerate regime. 


## Citation

The results obtained with this code are presented in the following paper:

> S. S. Sirletti, T. Namikawa, and B. D. Sherwin,
> *Alleviating CMB-BAO tension with a large cosmic birefringence angle during
> reionization induced by energy injection to axion-like particles*,
> [arXiv:2610.29764](https://arxiv.org/abs/2610.29764).

If you use this code in your work, we would be happy if you can cite:

```bibtex
@article{Sirletti:2026pgv,
    author = "Sirletti, Salvatore Samuele and Namikawa, Toshiya and Sherwin, Blake D.",
    title = "{Alleviating CMB-BAO tension with a large cosmic birefringence angle during reionization induced by energy injection to axion-like particles}",
    eprint = "2610.29764",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    month = "10",
    year = "2026"
}
```