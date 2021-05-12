#!/bin/bash
#
#SBATCH -J MTP
#SBATCH -A LiU-2019-26
#SBATCH -t 12:00:00
#SBATCH -N 2
#SBATCH --exclusive
#SBATCH -n 64
#
# Preamble
# consider putting mlp into repository
MLP_EXE=../../mlip-2/bin/mlp
TMP_DIR=./mtps_out
mkdir -p $TMP_DIR
NUM_POT=06
ELEMENT=Si
#NUM_TIMESTEPS=2000

# Body
# CHECK mtps_out and cfg_train, WAS USED FOR LOG SCALE
# This converts OUTCAR to the internal format .cfg
#i=10000
for i in 100
do
    $MLP_EXE train untrained_mtps/${NUM_POT}.mtp cfg_train/${ELEMENT}_train_${i}.cfg --trained-pot-name=$TMP_DIR/${ELEMENT}_${NUM_POT}_pot_${i}.mtp --max-iter=100
    #$MLP_EXE train untrained_mtps/${NUM_POT}.mtp cfg_train/${ELEMENT}_train_${i}.cfg --trained-pot-name=$TMP_DIR/${i}/${ELEMENT}_${NUM_POT}_pot_${i}_iter_${j}.mtp --max-iter=100
done
