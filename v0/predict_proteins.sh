#!/bin/bash

#SBATCH --job-name=prodigal
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=/gxfs_work1/geomar/smomw539/bird_bags/YuChens_pipeline//reports/prodigal_out
#SBATCH --error=/gxfs_work1/geomar/smomw539/bird_bags/YuChens_pipeline//reports/prodigal_err
#SBATCH --partition=cluster
#SBATCH --mail-type=END,FAIL,TIME-LIMIT
#SBATCH --mail-user=theokrueger.infomail@gmail.com

source $HOME/miniconda3/bin/activate prodigal

mkdir /gxfs_work1/geomar/smomw539/bird_bags/YuChens_pipeline/02_predicted/v0/
for bin in $(ls /gxfs_work1/geomar/smomw539/bird_bags/YuChens_pipeline/01_DNA_sequences/v0/); do
    prodigal \
        -i /gxfs_work1/geomar/smomw539/bird_bags/YuChens_pipeline/01_DNA_sequences/v0/$bin \
        -a /gxfs_work1/geomar/smomw539/bird_bags/YuChens_pipeline/02_predicted/v0/$bin.faa \
        -p meta
done
