#/bin/bash

# Preamble
# consider putting mlp into repository
MLP_EXE=../../mlip-2/bin/mlp
TMP_DIR=./mtps_out
mkdir -p $TMP_DIR
NUM_POT=06
ELEMENT=Al
#NUM_TIMESTEPS=2000

# Body
# CHECK mtps_out and cfg_train, WAS USED FOR LOG SCALE
# This converts OUTCAR to the internal format .cfg
for i in 180 190 200 210 220
do
    $MLP_EXE train untrained_mtps/${NUM_POT}.mtp cfg_train/${ELEMENT}_train_${i}.cfg --trained-pot-name=$TMP_DIR/${ELEMENT}_${NUM_POT}_pot_${i}.mtp --max-iter=100
done
