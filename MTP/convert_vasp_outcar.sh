#/bin/bash

# Preamble
# consider putting mlp into repository
MLP_EXE=../../mlip-2/bin/mlp
TMP_DIR=./cfg_out
mkdir -p $TMP_DIR

# Body
# This converts OUTCAR to the internal format .cfg
$MLP_EXE convert-cfg --input-format=vasp-outcar ../Si_300K/OUTCAR cfg_out/Si_relax.cfg
