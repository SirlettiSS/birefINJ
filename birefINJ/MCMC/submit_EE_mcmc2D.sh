#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH -A INF26_indark
#SBATCH --mem=64000
#SBATCH -p boost_usr_prod
#SBATCH --time 14:30:00
#SBATCH -J EE_mcmc2D
#SBATCH -o EE_mcmc2D.log
#SBATCH -e EE_mcmc2D.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=srlsvt@unife.it

cd /leonardo_work/INF26_indark/ssirlett/biref_from_injection/biref_iso_inject/

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

source ~/.bashrc
conda activate jupy_env

python EE_mcmc2D.py