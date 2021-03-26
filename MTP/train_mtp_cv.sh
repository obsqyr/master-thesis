#/bin/bash

# Preamble
# consider putting mlp into repository
MLP_EXE=../../mlip-2/bin/mlp
TMP_DIR=./mtps_out
mkdir -p $TMP_DIR
NUM_POT=06
ELEMENT=Al
NUM_TIMESTEPS=10

# Body
# This converts OUTCAR to the internal format .cfg
for i in 1 2 3 4 5 6 7 8 9 10
do
    $MLP_EXE train untrained_mtps/${NUM_POT}.mtp cfg_cv/train/${NUM_TIMESTEPS}/${ELEMENT}_traincv_${i}.cfg --trained-pot-name=cfg_cv/mtps_trained/${NUM_TIMESTEPS}/${ELEMENT}_${NUM_POT}_pot_${i}.mtp --max-iter=100
done
