#/bin/bash

# Preamble
# consider putting mlp into repository
MLP_EXE=../../mlip-2/bin/mlp
TMP_DIR=./mtps_out
mkdir -p $TMP_DIR
NUM_POT=14
ELEMENT=Al
#NUM_TIMESTEPS=2000

# Body
# This converts OUTCAR to the internal format .cfg
for i in 1000 2000 3000 4000 5000 6000 7000 8000 9000
do
    $MLP_EXE train untrained_mtps/${NUM_POT}.mtp cfg_train/${ELEMENT}_train_${i}.cfg --trained-pot-name=$TMP_DIR/${ELEMENT}_${NUM_POT}_pot_${i}.mtp --max-iter=100
done
