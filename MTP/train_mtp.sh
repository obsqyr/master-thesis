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
NUM_POT=10
ELEMENT=Al
#NUM_TIMESTEPS=2000

# Body
# CHECK mtps_out and cfg_train, WAS USED FOR LOG SCALE
# This converts OUTCAR to the internal format .cfg
#i = 100
for i in 100  
do
    for j in 0 1 2 3 4 5 6 7 8 9 
    do
    #$MLP_EXE train untrained_mtps/${NUM_POT}.mtp cfg_train/${ELEMENT}_train_${i}.cfg --trained-pot-name=$TMP_DIR/${ELEMENT}_${NUM_POT}_pot_${i}.mtp --max-iter=1000 --valid-cfgs=cfg_test/${ELEMENT}_test_1000.cfg | tee training_output/${ELEMENT}_${NUM_POT}_${i}.txt
	$MLP_EXE train untrained_mtps/${NUM_POT}.mtp cfg_train/${ELEMENT}_train_${i}.cfg --trained-pot-name=$TMP_DIR/${i}/${ELEMENT}_${NUM_POT}_pot_${i}_iter_${j}.mtp --max-iter=1000 --valid-cfgs=cfg_test/${ELEMENT}_test_1000.cfg | tee training_output_iter/${ELEMENT}_${NUM_POT}_${i}_${j}.txt
    done
done
