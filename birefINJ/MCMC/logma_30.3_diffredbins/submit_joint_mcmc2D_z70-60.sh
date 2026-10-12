#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH -A INF26_indark
#SBATCH --mem=64000
#SBATCH -p boost_usr_prod
#SBATCH --time 14:30:00
#SBATCH -J joint_z70-60
#SBATCH -o /leonardo_work/INF26_indark/ssirlett/biref_from_injection/biref_iso_inject/logma_30.3_diffredbins/joint_z70-60.log
#SBATCH -e /leonardo_work/INF26_indark/ssirlett/biref_from_injection/biref_iso_inject/logma_30.3_diffredbins/joint_z70-60.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=srlsvt@unife.it

cd /leonardo_work/INF26_indark/ssirlett/biref_from_injection/biref_iso_inject/logma_30.3_diffredbins/

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

export Z_START=70
export Z_END=60
export XI_MAX=8.89969e-03

source ~/.bashrc
conda activate jupy_env

python joint_mcmc2D.py