#!/bin/bash
#
#SBATCH -J testjob
#SBATCH -A LiU-2019-26
#SBATCH -t 48:00:00
#SBATCH -N 4
#SBATCH --exclusive
#SBATCH -n 128
#
./train_mtp.sh
