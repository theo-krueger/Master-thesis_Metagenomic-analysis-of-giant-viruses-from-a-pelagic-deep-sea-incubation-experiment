#!/bin/bash

#SBATCH --job-name=hmmscan_NCVOG
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=/gxfs_work1/geomar/smomw539/bird_bags/YuChens_pipeline/reports/hmmscan_NCVOG_out
#SBATCH --error=/gxfs_work1/geomar/smomw539/bird_bags/YuChens_pipeline/reports/hmmscan_NCVOG_err
#SBATCH --partition=cluster
#SBATCH --mail-type=END,FAIL,TIME-LIMIT
#SBATCH --mail-user=theokrueger.infomail@gmail.com

source $HOME/miniconda3/bin/activate hmmer

mkdir /gxfs_work1/geomar/smomw539/bird_bags/YuChens_pipeline/03_annotated/v0/
mkdir /gxfs_work1/geomar/smomw539/bird_bags/YuChens_pipeline/03_annotated/v0/NCVOGs

for bin in $(ls $base_dir/02_predicted/v2/); do
    hmmscan \
        -E 1e-04 \
        --tblout /gxfs_work1/geomar/smomw539/bird_bags/YuChens_pipeline/03_annotated/v0/NCVOGs/$bin.NCVOG_tblout \
        $WORK/bird_bags/db/Needham_2019/47_NCVOGs/47.proteins.megavirales_Needham_2019.hmms \
        /gxfs_work1/geomar/smomw539/bird_bags/YuChens_pipeline/02_predicted/v0/$bin &>/dev/null
done
