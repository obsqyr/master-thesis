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
# This converts OUTCAR to the internal format .cfg
for i in 100 200 300 400 500 600 700 800 900
do
    $MLP_EXE calc-errors mtps_out/${ELEMENT}_${NUM_POT}_pot_${i}.mtp cfg_test/${ELEMENT}_test_1000.cfg > test_results/${NUM_POT}/${ELEMENT}_${i}.txt
done
