#!/bin/bash
#
#SBATCH -J MTP
#SBATCH -A LiU-2019-26
#SBATCH -t 12:00:00
#SBATCH -N 2
#SBATCH --exclusive
#SBATCH -n 64
#
module load Python/3.8.3-anaconda-2020.07-extras-nsc1
source ~/tqtf33/bin/activate
python3 main_06_eq_0.py

echo "job completed"