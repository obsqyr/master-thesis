#!/bin/bash
#
#SBATCH -J testjob
#SBATCH -A LiU-2019-26
#SBATCH -t 48:00:00
#SBATCH -N 2
#SBATCH --exclusive
#SBATCH -n 64
#
./train_mtp.sh
