# The correct compilation of CLASS depends on the architecture used, e.g. it will be different if compiled 
# on a Macbook or on a HPC cluster. 

# For this motivation, here's a small code to compule birefclass on Leonardo HPC cluster. 

#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH -A INF26_indark
#SBATCH -p boost_usr_prod
#SBATCH --qos boost_qos_dbg
#SBATCH --time 00:30:00
#SBATCH -J recompile
#SBATCH -o recompile.log
#SBATCH -e recompile.err

module load gcc/12.2.0

cd /leonardo_work/INF26_indark/ssirlett/biref_from_injection/birefclass

# Restore new flags clean without ASAN
sed -i 's/OPTFLAG = .*/OPTFLAG = -O2 -fwrapv/g' Makefile
sed -i 's/CCFLAG = .*/CCFLAG = -g -fPIC/g' Makefile

# Recompile from the scratch
rm -rf build/ libclass.a class
make -j4 class

echo "Binary: $(ls -la class)"
echo "Testing CLASS:"
./class /leonardo_work/INF26_indark/ssirlett/biref_from_injection/pyrefCLASS/initial_file/biref_default_inject.ini ./cl_ref_biref.pre
echo "Exit code: $?"