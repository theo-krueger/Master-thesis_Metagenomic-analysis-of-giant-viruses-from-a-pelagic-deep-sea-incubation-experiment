#!/bin/bash

# predicting proteins
for bin in $(ls 01_DNA_sequences/); do
sbatch --job-name=$bin.prodigal \
--nodes=1 \
--tasks-per-node=32 \
--cpus-per-task=1 \
--mem=120G \
--time=01:00:00 \
--output=reports/prodigal_out \
--error=reports/prodigal_err \
--wrap="\
prodigal -i 01_DNA_sequences/$bin -a 02_predicted/$bin.faa -p meta"
done