#/bin/bash

# Preamble
# consider putting mlp into repository
MLP_EXE=../../mlip-2/bin/mlp
TMP_DIR=./mtps_out_eq
mkdir -p $TMP_DIR
NUM_POT=06
ELEMENT=Al
#NUM_TIMESTEPS=2000
EQ=2000

# Body
# CHECK mtps_out and test_results, WAS USED FOR LOG SCALE
# This converts OUTCAR to the internal format .cfg
for i in 200
do
    $MLP_EXE calc-errors mtps_out_eq/${ELEMENT}_${NUM_POT}_pot_${i}_eq_${EQ}.mtp cfg_test/${ELEMENT}_test_1000.cfg > test_results_eq/${NUM_POT}/${ELEMENT}_${i}_eq_${EQ}.txt
done
