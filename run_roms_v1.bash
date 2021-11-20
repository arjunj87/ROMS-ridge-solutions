#!/bin/bash
#SBATCH -A cla119
#SBATCH --job-name="ellipse_3D_infinite"  
#SBATCH --output="ellipse_3D_infinite.%j.%N.out"  
#SBATCH --partition=compute  
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=24  
#SBATCH --export=ALL  
#SBATCH -t 16:30:00  

# Start of user commands - comments start with a hash sign (#)

cd /home/arjunj87/roms_ellipticalridge/Racetrack_Bu_3dot2/ 
ibrun -np 192 ./roms island.in 
