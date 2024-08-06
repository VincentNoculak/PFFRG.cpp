#!/bin/bash

#SBATCH --job-name=Triangular      # Job name, will show up in squeue Out-put
#SBATCH --ntasks=128                     # Number of cores
#SBATCH --nodes=1                      # Ensure that all cores are on one machine
#SBATCH --time=1-00:00:00                # Runtime in DAYS-HH:MM:SS format
#SBATCH --mem-per-cpu=600             # Memory per cpu in MB (see also --mem) 
#SBATCH --output=PKo.dat   # File to which standard Out- will be written
#SBATCH --error=PKe.dat    # File to which standard err will be written
#SBATCH --mail-type=END                # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=v.noculak@fu-berlin.de   # Email to which notifications will be sent 
#SBATCH -p normal

export OMP_NUM_THREADS=128

# run your program...
./PFFRG "$1"

