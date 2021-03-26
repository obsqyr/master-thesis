#/bin/bash

# Preamble
# consider putting mlp into repository
MLP_EXE=../../mlip-2/bin/mlp
TMP_DIR=./mtps_out_log
mkdir -p $TMP_DIR
NUM_POT=06
ELEMENT=Si
#NUM_TIMESTEPS=2000

# Body
# CHECK mtps_out and test_results, WAS USED FOR LOG SCALE
# This converts OUTCAR to the internal format .cfg
for i in 2 3 4 5 6 7 8 9 10 11 13 14 16 17 19 22 24 27 30 34 38 42 47 53 59 66 73 82 92 103 115 129 144 161 180 202 226 253 283 317
do
    $MLP_EXE calc-errors mtps_out_log/${ELEMENT}_${NUM_POT}_pot_${i}.mtp cfg_train_log/${ELEMENT}_train_${i}.cfg > train_results_log/${NUM_POT}/${ELEMENT}_${i}.txt
done
