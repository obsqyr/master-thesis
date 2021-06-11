#/bin/bash

# Preamble
# consider putting mlp into repository
MLP_EXE=../../mlip-2/bin/mlp
TMP_DIR=./mtps_out_eq
mkdir -p $TMP_DIR
NUM_POT=10
ELEMENT=Si
EQ=2000
#NUM_TIMESTEPS=2000

# Body
# CHECK mtps_out and cfg_train, WAS USED FOR LOG SCALE
# This converts OUTCAR to the internal format .cfg
#i=10000
for i in 1 10 100 1000 10000
do
    $MLP_EXE train untrained_mtps/${NUM_POT}.mtp cfg_train_eq/${ELEMENT}_train_${i}_eq_${EQ}.cfg --trained-pot-name=$TMP_DIR/${ELEMENT}_${NUM_POT}_pot_${i}_eq_${EQ}.mtp --max-iter=1000 | tee training_output/${ELEMENT}_${NUM_POT}_eq_${EQ}_${i}.txt
    #$MLP_EXE train untrained_mtps/${NUM_POT}.mtp cfg_train_eq/${ELEMENT}_train_${i}_eq_${EQ}.cfg --trained-pot-name=$TMP_DIR/${i}/${ELEMENT}_${NUM_POT}_pot_${i}_eq_${EQ}_iter_${j}.mtp --max-iter=100
done
