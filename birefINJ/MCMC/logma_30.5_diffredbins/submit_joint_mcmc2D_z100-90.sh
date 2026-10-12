#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH -A INF26_indark
#SBATCH --mem=64000
#SBATCH -p boost_usr_prod
#SBATCH --time 14:30:00
#SBATCH -J joint_z100-90
#SBATCH -o /leonardo_work/INF26_indark/ssirlett/biref_from_injection/biref_iso_inject/logma_30.5_diffredbins/joint_z100-90.log
#SBATCH -e /leonardo_work/INF26_indark/ssirlett/biref_from_injection/biref_iso_inject/logma_30.5_diffredbins/joint_z100-90.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=srlsvt@unife.it

cd /leonardo_work/INF26_indark/ssirlett/biref_from_injection/biref_iso_inject/logma_30.5_diffredbins/

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

export Z_START=100
export Z_END=90
export XI_MAX=2.01712e-02

source ~/.bashrc
conda activate jupy_env

python joint_mcmc2D.py