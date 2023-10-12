#!/bin/sh

#SBATCH --job-name=diamond_pdb
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --output=reports/diamond_pdb_out
#SBATCH --error=reports/diamond_pdb_err
#SBATCH --partition=cluster
#SBATCH --mail-type=END,FAIL,TIME-LIMIT
#SBATCH --mail-user=theokrueger.infomail@gmail.com

source $HOME/miniconda3/bin/activate diamond-2.1.4

pred_dir=$WORK/bird_bags/YuChens_pipeline/02_predicted/v2_ESOM
out_dir=$WORK/bird_bags/YuChens_pipeline/03_annotated/v2_ESOM

for file in $(ls $pred_dir); do
    diamond blastp \
        --query $pred_dir/$file \
        --db $WORK/bird_bags/db/pdb_230517/pdbaa.gz \
        --out $out_dir/pdb/$file.vs.pdb.diamond \
        --threads 32 \
        --evalue 1e-04 \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle staxids sscinames \
        --max-target-seqs 10 
done
