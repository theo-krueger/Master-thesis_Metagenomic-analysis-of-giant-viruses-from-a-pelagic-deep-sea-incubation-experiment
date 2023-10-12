#!/bin/bash

# predicting proteins
for bin in $(ls 01_DNA_sequences/v2/); do
sbatch --job-name=$bin.prodigal \
--nodes=1 \
--tasks-per-node=32 \
--cpus-per-task=1 \
--mem=120G \
--time=01:00:00 \
--output=reports/prodigal_out \
--error=reports/prodigal_err \
--partition=cluster \
--mail-type=END,FAIL,TIME-LIMIT \
--wrap="\
source $HOME/miniconda3/bin/activate prodigal
prodigal -i 01_DNA_sequences/v2/$bin -a 02_predicted/v2.1v2.1/$bin.faa -p meta"
done

# one file got lost somehow
sbatch --job-name=coAssem_concoct_77_4.prodigal \
--nodes=1 \
--tasks-per-node=32 \
--cpus-per-task=1 \
--mem=120G \
--time=01:00:00 \
--output=reports/prodigal_77_4_out \
--error=reports/prodigal_77_4_err \
--partition=cluster \
--mail-type=END,FAIL,TIME-LIMIT \
--wrap="\
source $HOME/miniconda3/bin/activate prodigal
prodigal -i 01_DNA_sequences/v2/BB_coAssem_output_concoct_77_4-contigs.fa -a 02_predicted/v2/BB_coAssem_output_concoct_77_4-contigs.fa.faa -p meta"

# bins that went down with time
for bin in $(ls $WORK/bird_bags/YuChens_pipeline/01_DNA_sequences/v2.2/); do
sbatch \
    --job-name=$bin.prodigal \
    --nodes=1 \
    --tasks-per-node=32 \
    --cpus-per-task=1 \
    --mem=120G \
    --time=01:00:00 \
    --output=$WORK/bird_bags/YuChens_pipeline/reports/prodigal_out \
    --error=$WORK/bird_bags/YuChens_pipeline/reports/prodigal_err \
    --partition=cluster \
    --mail-type=END,FAIL,TIME-LIMIT \
    --wrap="\
    source $HOME/miniconda3/bin/activate prodigal
    prodigal -i $WORK/bird_bags/YuChens_pipeline/01_DNA_sequences/v2.2/$bin -a $WORK/bird_bags/YuChens_pipeline/02_predicted/v2.2/$bin.faa -p meta
    "
done

# for rest of polB bins
conda activate prodigal
for bin in $(ls $WORK/bird_bags/YuChens_pipeline/01_DNA_sequences/v2_polB_bins/); do
    prodigal \
        -i $WORK/bird_bags/YuChens_pipeline/01_DNA_sequences/v2_polB_bins/$bin \
        -a $WORK/bird_bags/YuChens_pipeline/02_predicted/v2_polB_bins/$bin.faa \
        -p meta
done
