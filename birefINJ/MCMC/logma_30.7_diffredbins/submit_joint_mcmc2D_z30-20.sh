#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH -A INF26_indark
#SBATCH --mem=64000
#SBATCH -p boost_usr_prod
#SBATCH --time 14:30:00
#SBATCH -J joint_z30-20
#SBATCH -o /leonardo_work/INF26_indark/ssirlett/biref_from_injection/biref_iso_inject/logma_30.7_diffredbins/joint_z30-20.log
#SBATCH -e /leonardo_work/INF26_indark/ssirlett/biref_from_injection/biref_iso_inject/logma_30.7_diffredbins/joint_z30-20.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=srlsvt@unife.it

cd /leonardo_work/INF26_indark/ssirlett/biref_from_injection/biref_iso_inject/logma_30.7_diffredbins/

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

export Z_START=30
export Z_END=20
export XI_MAX=2.77233e-02

source ~/.bashrc
conda activate jupy_env

python joint_mcmc2D.py