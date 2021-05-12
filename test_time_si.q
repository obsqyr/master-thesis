#!/bin/bash
#
#SBATCH -J MTP
#SBATCH -A LiU-2019-26
#SBATCH -t 12:00:00
#SBATCH --exclusive
#SBATCH -n 1
#
module load Python/3.8.3-anaconda-2020.07-extras-nsc1
source ~/tqtf33/bin/activate
time python3 test_time_si.py

echo "job completed"
