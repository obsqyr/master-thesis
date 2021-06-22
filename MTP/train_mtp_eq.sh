#/bin/bash

# Preamble
# consider putting mlp into repository
MLP_EXE=../../mlip-2/bin/mlp
TMP_DIR=./mtps_out_eq
mkdir -p $TMP_DIR
NUM_POT=06
ELEMENT=Al
EQ=2000
#NUM_TIMESTEPS=2000

# Body
# CHECK mtps_out and cfg_train, WAS USED FOR LOG SCALE
# This converts OUTCAR to the internal format .cfg
i=100
for j in 0 1 2 3 4 5 6 7 8 9
do
    #$MLP_EXE train untrained_mtps/${NUM_POT}.mtp cfg_train_eq/${ELEMENT}_train_${i}_eq_${EQ}.cfg --trained-pot-name=$TMP_DIR/${ELEMENT}_${NUM_POT}_pot_${i}_eq_${EQ}.mtp --max-iter=1000 | tee training_output/${ELEMENT}_${NUM_POT}_eq_${EQ}_${i}.txt
    $MLP_EXE train untrained_mtps/${NUM_POT}.mtp cfg_train_eq/${ELEMENT}_train_${i}_eq_${EQ}.cfg --trained-pot-name=$TMP_DIR/${i}/${ELEMENT}_${NUM_POT}_pot_${i}_eq_${EQ}_iter_${j}.mtp --max-iter=1000 --valid-cfgs=cfg_test/${ELEMENT}_test_1000.cfg | tee training_output_iter/${ELEMENT}_${NUM_POT}_${i}_${j}.txt
done
