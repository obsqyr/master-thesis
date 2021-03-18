#!/bin/sh

# Preamble
# consider putting mlp into repository
MLP_EXE=../mlip-2/bin/mlp
TMP_DIR=./mtps_out
mkdir -p $TMP_DIR
NUM_POT=06
ELEMENT=Al
i=1000

# Body
# $1: path to mtp
$MLP_EXE calc-efs ${1} MTP/atom.cfg MTP/pred.cfg

#$MLP_EXE calc-efs MTP/mtps_out/${ELEMENT}_${NUM_POT}_pot_${i}.mtp MTP/atom.cfg MTP/pred.cfg

#cfg_predicted/${ELEMENT}_${NUM_POT}_1.cfg 
#$MLP_EXE calc-efs mtps_out/Al_06_pot_1000.mtp atom.cfg pred_test.cfg
