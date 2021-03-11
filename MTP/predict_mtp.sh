#/bin/bash

# Preamble
# consider putting mlp into repository
MLP_EXE=../../mlip-2/bin/mlp
TMP_DIR=./mtps_out
mkdir -p $TMP_DIR
NUM_POT=06
ELEMENT=Al
i=1000

# Body
# This converts OUTCAR to the internal format .cfg
$MLP_EXE calc-efs mtps_out/${ELEMENT}_${NUM_POT}_pot_${i}.mtp atom.cfg cfg_predicted/${ELEMENT}_${NUM_POT}_1.cfg 
#$MLP_EXE calc-efs mtps_out/Al_06_pot_1000.mtp cfg_train/Al_train_1000.cfg pred_test.cfg
