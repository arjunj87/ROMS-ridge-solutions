#!/bin/bash
#SBATCH -A cla119 
#SBATCH --job-name="seamount"  
#SBATCH --output="seamount.%j.%N.out"  
#SBATCH --partition=compute  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24  
#SBATCH --export=ALL  
#SBATCH -t 00:40:00  

# Start of user commands - comments start with a hash sign (#)

./ncjoin elliptical_his.****.???.nc
./ncjoin elliptical_diags_uv.0000.???.nc
