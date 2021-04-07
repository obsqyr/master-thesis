#/bin/bash

# Preamble
# consider putting mlp into repository
MLP_EXE=../../mlip-2/bin/mlp
TMP_DIR=./mtps_out
mkdir -p $TMP_DIR
NUM_POT=10
ELEMENT=Si
#NUM_TIMESTEPS=2000

# Body
# CHECK mtps_out and test_results, WAS USED FOR LOG SCALE
# This converts OUTCAR to the internal format .cfg
for i in 180 190 200 210 220 
do
    $MLP_EXE calc-errors mtps_out/${ELEMENT}_${NUM_POT}_pot_${i}.mtp cfg_train/${ELEMENT}_train_${i}.cfg > test_results/${NUM_POT}/${ELEMENT}_${i}.txt
done
